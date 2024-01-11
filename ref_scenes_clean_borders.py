"""
Apply border masking to reference scene Assets.

Andrew Tedstone, July 2023
"""

# Generic libraries
import argparse
import tomli
import re

# Earth Engine
import ee
import geemap
ee.Initialize()

from load_env import *
from gee_functions import * 
from gee_s1 import *

debug = False

if __name__ == '__main__':    

    RETR_MAX = 500

    p = argparse.ArgumentParser(description='Apply border masking to reference scene assets.')
    p.add_argument('config_file', type=str, help='str, path to job configuration TOML file')
    p.add_argument('-debug', dest='debug', default=False, action='store_true', help='Run in debug mode')
    args = p.parse_args()

    with open(args.config_file, "rb") as f:
        conf = tomli.load(f)

    if conf['dfilter'] is not None:
        despeckle_suffix = '_{f}{p}'.format(f=conf['dfilter'], p=conf['filterpx'])
    match_pattern = r'ref_{m}_{p}_[0-9]+_{ys}_{ye}_{ds}_{de}_{res}m{despeck}_WITHOUTBORDERCLEAN$'.format(
        m=conf['mode'], 
        p=conf['pol'],
        ds=conf['ref_start_doy'],
        de=conf['ref_end_doy'], 
        ys=conf['ref_start_year'],
        ye=conf['ref_end_year'],
        res=conf['export_resolution_m'],
        despeck=despeckle_suffix)
    ref_assets = ee.data.listAssets({'parent':'projects/firn-retention/assets/'})
    uris = []
    for asset in ref_assets['assets']:
        if re.match(match_pattern, asset['id'].split('/')[-1]) is not None:
            print(asset['id'])
            uris.append(asset['id'])

    def load_im(uri):
        i = ee.Image(uri)
        return i.updateMask(i.neq(-9999))
    refs = ee.ImageCollection([load_im(i) for i in uris])

    def clean_borders(im):
        """
        # Find edges to constrain to. This approach is based on: 
        # https://gis.stackexchange.com/questions/347990/masking-edges-of-s1-composite-in-google-earth-engine
        # but is modified to work with a median image rather than a stack of raw s1 scenes.
        # Distance to nearest non-zero value (within 128px neighbourhood).

        THIS ASSUMES 100PX ASSETS!!
        """
        entropy = im.select(conf['pol']).toInt().entropy(ee.Kernel.rectangle(3,3, units='pixels'))
        pixelsToMask = im.mask().Not().fastDistanceTransform(128, 'pixels').sqrt()
        metersToMask = pixelsToMask.multiply(ee.Image.pixelArea().sqrt()).rename('metersToMask')
        notBorder = metersToMask.gte(500).And(pixelsToMask.gt(2))
        border = metersToMask.lte(500).And(pixelsToMask.lte(3))
        bad_border = border.And(entropy.gt(1.0))
        # Deal with small artefacts not connected to the main image area
        connected = bad_border.Not().connectedPixelCount(30, eightConnected=False)
        im_masked = im.updateMask(bad_border.Not().And(connected.gte(30)))
        return im_masked

    refs = refs.map(clean_borders)


    def export_ref(im, region):
        kws = dict(crs=conf['crs'], scale=conf['export_resolution_m'], maxPixels=1e13)
        fn_prefix = 'ref_{mode}_{pol}_{ron}_{ys}_{ye}_{ds}_{de}_{res}m{d}'
        im_info = im.getInfo()
        fn_prefix = fn_prefix.format(
            mode=im_info['properties']['mode'],
            pol=im_info['bands'][0]['id'],
            ron=im_info['properties']['orbit'],
            ys=im_info['properties']['year_start'], 
            ye=im_info['properties']['year_end'], 
            ds=conf['ref_start_doy'],
            de=conf['ref_end_doy'],
            res=conf['export_resolution_m'],
            d=despeckle_suffix)
        task = ee.batch.Export.image.toAsset(
            image=im.unmask(-9999),
            region=region,
            description=fn_prefix,
            assetId='{a}/{p}'.format(a=conf['assets_path'], p=fn_prefix), 
            **kws)
        task.start()
        return task

    print('Calculating number of scenes to export...')
    nrefs = refs.size().getInfo()
    
    offset = 0
    print('{n} reference scenes to export ...'.format(n=nrefs))
    

    img_lst = refs.toList(RETR_MAX, offset=offset)
    imgs_in_lst = img_lst.size().getInfo()
    task_store = []
    count = 0

    for i in range(0, imgs_in_lst):
        image = ee.Image(img_lst.get(i))

        t = export_ref(image.select(conf['pol']), image.geometry())

        status = t.status()
        print('{n}/{t} {d} : {s}'.format(n=count+1, t=imgs_in_lst, 
                d=status['description'], s=status['state']))

        task_store.append(t)
        count += 1


