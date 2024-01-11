"""
Sentinel-1 GEE functions

AT September 2022
"""
import ee
import sys
sys.path.append("/Users/tedstona/scripts/gee_s1_ard/python-api")
import speckle_filter

## Sentinel-1
NOISE_FLOOR_DB = -25

NATIVE_RES = {
    'IW': 25,
    'EW': 40
}

def toDb(linear, bands):
    """
    https://gis.stackexchange.com/questions/424225/convert-sentinel-1-images-data-from-db-to-linear
    """
    return linear.addBands(
    ee.Image().expression('10 * log10(linear)', {'linear':linear.select(bands)}).rename(bands),
    overwrite=True # Replace the bands to keep image properties
    )


def import_collection(
    coll, 
    mode, 
    pol, 
    do_remove_dark_borders=False,
    do_mask_floor=False
    ):
    """
    Import and filter an ImageCollection of Sentinel-1 acquisitions.
    """
    # Get Sentinel 1 acquisitions of interest
    s1 = ee.ImageCollection(coll) \
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', pol)) 

    s1 = s1.filter(ee.Filter.eq('instrumentMode', mode))
    if do_remove_dark_borders:
        # This has to happen before any other operation, because it assumes a homogeneous (filled) mask
        s1 = s1.map(mask_border) #remove_dark_borders
    if do_mask_floor:
        s1 = s1.map(remove_below_noise_floor)

    return s1

def force_native_projection(im, mode):
    """
    im
    mode : IW, EW, etc (see NATIVE_RES)

    Use by writing a wrapping function, e.g.:

    def force_native_projection_iw(im):
        return force_native_projection(im, mode='IW')
    my_collection.map(force_native_projection_iw)

    """
    return im.reproject(im.projection(), scale=NATIVE_RES[mode])


def prepare_despeckle(fname, fpx, s1_mode, band,
    force_reproject=False):
    """
    Despeckling needs to be done at the native resolution of the imaging mode!
    """

    # def _project_im(im):
    #     return im.reproject(im.select(band).projection(), scale=NATIVE_RES[s1_mode])
    
    def despeckle_box(im):
        """
        This is in PIXELS, so sensible selection of this depends on the export resolution.
        Liang used 9*9 boxcar at 40 m resolution = 360 m.
        So gridding here at 100 m, a 3 px boxcar makes sense
        """
        return speckle_filter.boxcar(im, fpx)
    
    def despeckle_lee(im):
        return speckle_filter.leefilter(im, fpx) #_project_im(im)

    if fname == 'lee':
        return despeckle_lee
    elif fname == 'box':
        return despeckle_box
    else:
        raise ValueError('The filter you have specified is not implemented.')



def ang_corr(im):
    """ Following eq. 2 of Johnson et al., 2020, Remote Sensing of Env. """
    i = im.expression("b('HH') - (b('angle') - 32) * -0.17")
    i = i.rename('HH_corr')
    im = im.addBands(srcImg=i)
    return im


def get_orbit(im):
    """ 
    Get relativeOrbitNumber of an Image. Returns Feature.

    Have to return computed objects as Features or Images.
    Features have to have a geometry, here let's spoof one.
    """
    try:
        return ee.Feature(ee.Geometry.Point(0,0), {
            'orbit':im.get('relativeOrbitNumber_start'), 
            'mode':im.get('instrumentMode')
            })
    except:
        return ee.Feature(ee.Geometry.Point(0,0), {
            'orbit':0
            })

    
def create_annual_ref_scenes(collection, mode,
    year_start, year_end,
    ref_doy_start=1, ref_doy_end=90):
    """  DEPRECATED !! See create_ref_scenes.

    Create reference mosaics from winter-time data, one per orbit for each year.
    
    Returns: ee.ImageCollection of reference scenes.
    Metadata added to return: year, orbit, nscenes.
    
    """

    # Filter the provided collection to only the reference period
    ref_collection = collection.filter(ee.Filter.dayOfYear(start=ref_doy_start, end=ref_doy_end)) \
        .filter(ee.Filter.eq('instrumentMode', mode))

    def _create_refs(year):
        """ Create reference images, one per orbit, for a specified year.
        
        This is a private function because it uses the filtered ImageCollection from the parent function scope.
        
        Returns: ee.ImageCollection
        """
        
        def _create_orbit_ref(orbit):
            # Reduce to an image for the orbit this year
            ims_in_orbit = ref_year_collection.filter(ee.Filter.eq('relativeOrbitNumber_start', orbit))
            ref_this_orbit = ims_in_orbit.mean()
            # Save some associated metadata - needed for filtering by user
            n_ims = ims_in_orbit.size()
            ref_this_orbit = ref_this_orbit.set('year', year).set('orbit', orbit).set('nscenes', n_ims).set('mode', mode)
            # Return an ee.Image
            return ref_this_orbit

        # Filter ImageCollection to given year
        ref_year_collection = ref_collection.filter(ee.Filter.calendarRange(year, year, 'year'))
        
        # Find all the orbits in the year's ImageCollection
        orbits = ref_year_collection.map(get_orbit)
        orbits = orbits.reduceColumns(reducer=ee.Reducer.toList(), selectors=['orbit'])        
        unique_orbits = ee.List.distinct(orbits.get('list'))
        
        refs = unique_orbits.map(_create_orbit_ref)
        refs = ee.ImageCollection.fromImages(refs)
        return refs
    
    years = ee.List.sequence(year_start, year_end)
    
    # https://gis.stackexchange.com/questions/423392/merge-a-list-of-imagecollection-into-a-single-imagecollection-in-google-earth-en/423461
    ref_collection = ee.ImageCollection( # Cast result to ImageCollection
        ee.FeatureCollection( # Create an intermediate FeatureCollection containing ImageCollections
            years.map(_create_refs)
        ).flatten() # Remove the intermediate feature collection with the merged image collections
    )
    
    return ref_collection


def create_ref_scenes(collection, mode, pol,
    year_start, year_end,
    ref_doy_start=1, ref_doy_end=90,
    blend_borders=False, blend_borders_nibble_px=4):
    """ Create reference mosaics from winter-time data, one per orbit, using all data in specified time range.
    
    Returns: ee.ImageCollection of reference scenes.
    Metadata added to return: year, orbit, nscenes.

    blend_borders_nibble_px: for 100 m export, 4 px works well. Equiv. to 400 m window.
    This suggests that for 20 m export, 400/20 = 20 pixels.
    
    """

    # Filter the provided collection to only the reference period
    ref_collection = collection.filter(ee.Filter.dayOfYear(start=ref_doy_start, end=ref_doy_end)) \
        .filter(ee.Filter.eq('instrumentMode', mode)) \
        .filter(ee.Filter.calendarRange(year_start, year_end, 'year'))

    if blend_borders:
        # We first save the default mask to the image.
        def save_orig_mask(im):
            m = im.select(pol).mask()
            m = m.rename('mask_orig')
            return im.addBands(m)
        ref_collection = ref_collection.map(save_orig_mask)

        # Next we nibble the valid image area back,  #####in testing 4 px works well if exporting at 100 m res.
        def nibble_px(im):
            return nibble_borders(im, r=blend_borders_nibble_px)
        ref_collection = ref_collection.map(nibble_px)

        # Then we find the pixels that we want to fill later 
        # While create_ref_scenes() does not apply the filling, it does need to compute max-based
        # masks that apply for the whole reference scene, hence mask_orig and difference need to be
        # computed/added to each image before create_ref_scenes() runs.
        def identify_fill(im):
            a = im.select('mask_orig').unmask().subtract(im.select(pol).mask())
            a = a.rename('difference')
            im = im.addBands(a)
            return im
        ref_collection = ref_collection.map(identify_fill)
        
    def _create_orbit_ref(orbit):
        # Reduce to an image for the orbit
        ims_in_orbit = ref_collection.filter(ee.Filter.eq('relativeOrbitNumber_start', orbit))

        if blend_borders:
            diffmask = ims_in_orbit.select('difference').max()
            diffmask = diffmask.rename('difference')
            origmask = ims_in_orbit.select('mask_orig').max()
            origmask = origmask.rename('mask_orig')

        ref_this_orbit = ims_in_orbit.select(pol).median()
        if blend_borders:
            ref_this_orbit = ref_this_orbit.addBands(diffmask, overwrite=True)
            ref_this_orbit = ref_this_orbit.addBands(origmask, overwrite=True)

        # Save some associated metadata - needed for filtering by user
        n_ims = ims_in_orbit.size()
        pass_type = ims_in_orbit.first().get('orbitProperties_pass')
        # This block below causes very slow exports...disabled
        # date_start = ims_in_orbit.sort('system:time_start', True).first().date().millis()
        # date_end = ims_in_orbit.sort('system:time_start',  False).first().date().millis()
        # .set('system:time_start', date_start)
        # .set('system:time_end', date_end)
        ref_this_orbit = ( ref_this_orbit
            .set('orbit', orbit)
            .set('nscenes', n_ims)
            .set('mode', mode)
            .set('year_start', year_start)
            .set('year_end', year_end)
            .set('orbitProperties_pass', pass_type)
        )
        # Return an ee.Image
        #return filled
        
        return ref_this_orbit
    
    # Find all the orbits in the year's ImageCollection
    orbits = ref_collection.map(get_orbit)
    orbits = orbits.reduceColumns(reducer=ee.Reducer.toList(), selectors=['orbit'])        
    unique_orbits = ee.List.distinct(orbits.get('list'))
    
    #unique_orbits = ee.List(list(range(1,176)))
    refs = unique_orbits.map(_create_orbit_ref)
    refs = ee.ImageCollection.fromImages(refs)

    if blend_borders:
        # Now we fill the areas of interesting using a focal median filter.
        def med_mask(im):
            # unmask the ref scene, re-apply the original mask and compute the area median
            fill = im.unmask().updateMask(im.select('mask_orig')).focalMedian(radius=10, units='pixels')
            # blend the filtered image with the 'proper' image.
            filled = im.blend(fill.mask(im.unmask().select('difference')))
            return filled
        refs = refs.map(med_mask)

    # These steps have to be applied at a sensible resolution! Let's use the native resolution.
    def _project_im(im):
        return im.reproject(im.select(pol).projection(), scale=NATIVE_RES[mode])
    refs = refs.map(_project_im)


    return refs
    
def _debias_collection(
    im, 
    ref, 
    band, 
    noise_floor=NOISE_FLOOR_DB,
    apply_ref_mask=True
    ):
    """ Debias an image using a reference image """
    result = im.select(band).subtract(ref.select(band))
    # Apply a minimum backscatter value threshold, following Johnson et al. 2020 Remote Sensing of Env.
    result = remove_below_noise_floor(result, noise_floor=noise_floor)
    # Apply the reference scene mask only after the reference has been used !! (2023-07-04)
    if apply_ref_mask:
        result = result.updateMask(ref.neq(-9999))
    # These have to come after the above lines, otherwise GEE mistakes that the result is an ee.Element.
    result = result.copyProperties(im, exclude=['reference'])
    result = result.copyProperties(im, properties=['system:time_start'])
    return result

def debias_collection_via_join(im, band, noise_floor=NOISE_FLOOR_DB):
    """ Debias an image which contains a sub-image named 'reference' following a join operation """
    summer = ee.Image(im)
    ref = ee.Image(summer.get('reference')) 
    return _debias_collection(summer, ref, band)
    

def debias_collection_via_lookup(im, band, refs_valid, noise_floor=NOISE_FLOOR_DB):
    """ 
    Implemented as work-around to GEE bug concerning application of join.
    Requires:
    - a refs_valid ImageCollection
    - a properly-filtered scenes ImageCollection (e.g. by orbits list)
    """
    summer = ee.Image(im)
    ref = refs_valid.filter(ee.Filter.eq('orbit', summer.get('relativeOrbitNumber_start'))).first()
    return _debias_collection(summer, ref, band)


def remove_below_noise_floor(im, noise_floor=NOISE_FLOOR_DB):
    """ Mask out data which fall below the specified noise floor. """
    mod = im.gt(noise_floor)
    return im.updateMask(mod)


def remove_dark_borders(im):
    """
    Probably supply only single-band images.

    Follows approach outlined in Hu et al. (2022), DOI: 10.1109/JSTARS.2022.3192409.
    """
    # First shrink the existing image mask by the specified distance
    msk_shrunk = im.mask().focalMin(radius=10000, units='meters', iterations=1, kernelType='square')
    # Then define the area of interest as the differnce between the shrunk mask and the original mask grown by 1 km
    aoi = im.mask() \
            .gt(msk_shrunk)
        #.focalMax(radius=1000, units='meters', iterations=1, kernelType='square') \
    # In this AOI, mask out values with a value less than the specified value.
    mask_values = (im.lt(-10).And(aoi.gt(0))).Not()
    # Apply mask to image and return in.
    return im.updateMask(mask_values)


def nibble_borders(im, r=2):
    m = (im.mask()).focalMin(radius=r, units='pixels', iterations=1)
    return im.updateMask(m)


def fill_blend(im, px=15):
    filled = im.focalMedian(px, 'square', 'pixels', 2)
    blended = filled.blend(im)
    return blended 

def mask_border(image, band='HV'):
    """
    https://code.earthengine.google.com/97698941ca6e15516ae890a0c2676e0c
    https://gis.stackexchange.com/questions/347990/masking-edges-of-s1-composite-in-google-earth-engine
    """
    totalSlices = ee.Number(image.get('totalSlices'))
    sliceNumber = ee.Number(image.get('sliceNumber'))
    middleSlice = ee.Image(sliceNumber.gt(1).And(sliceNumber.lt(totalSlices)))
    mask = image.select(band).mask().reduce(ee.Reducer.min()).floor()
    pixelsToMask = mask.Not().fastDistanceTransform(128, 'pixels').sqrt()
    metersToMask = pixelsToMask.multiply(ee.Image.pixelArea().sqrt()).rename('metersToMask')
    notBorder = metersToMask.gte(500).And(pixelsToMask.gt(2))
    angle = image.select('angle')
    return image.updateMask(
        angle.gt(31).And(angle.lt(45))
        .And(middleSlice.Or(notBorder))
        )    

    