"""
Debias summer scenes using winter references, threshold to identify melt,
then export end-of-season pan-GrIS metrics.

Andrew Tedstone, Sept 2022, based on scripts Nov 2021 to Aug 2022.
"""

# Generic libraries
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import os
import geopandas as gpd
import sys
import tomli
import re

# Earth Engine
# import ee
# import geemap
# ee.Initialize()
# sys.path.append("/Users/tedstona/scripts/gee_s1_ard/python-api")
# import speckle_filter

from load_env import *
from load_gee_env import *
from gee_functions import *
from gee_s1 import *

#############
debias_method = 'bylist' #valid options: join, bylist
#############

config_file_fn = sys.argv[1]

with open(config_file_fn, "rb") as f:
    conf = tomli.load(f)

if conf['ref_type'] == 'dynamic':
    use_ref_assets = False
elif conf['ref_type'] == 'asset':
    use_ref_assets = True
else:
    raise ValueError('Unknown ref_type setting.')

## Load scenes
s1_preproc = import_collection(conf['collection'], conf['mode'], conf['pol'], 
    do_remove_dark_borders=True,
    do_mask_floor=False)

## Load reference scenes with basic pattern matching - they are then filtered more below with metadata.
if use_ref_assets: #r'ref_{m}_{p}_[0-9]+_[0-9]+_[0-9]+_{ds}_{de}$'.format(
    if conf['dfilter'] is not None:
        despeck = '_{f}{p}'.format(f=conf['dfilter'], p=conf['filterpx'])
    match_pattern = r'ref_{m}_{p}_[0-9]+_{ys}_{ye}_{ds}_{de}_{res}m{despeck}$'.format(
        m=conf['mode'], 
        p=conf['pol'],
        ds=conf['ref_start_doy'],
        de=conf['ref_end_doy'], 
        ys=conf['ref_start_year'],
        ye=conf['ref_end_year'],
        res=conf['export_resolution_m'],
        despeck=despeck)
    ref_assets = ee.data.listAssets({'parent':'projects/firn-retention/assets/'})
    uris = []
    for asset in ref_assets['assets']:
        #if 'ref' in asset['id']:
        if re.match(match_pattern, asset['id'].split('/')[-1]) is not None:
            print(asset['id'])
            uris.append(asset['id'])

    def load_im(uri):
        i = ee.Image(uri)
        return i
        #return i.updateMask(i.neq(-9999))
    refs_all = ee.ImageCollection([load_im(i) for i in uris])
else:
    print('Creating ref scenes dynamically (WARNING: always multi-year!) (WARNING: no despeckling)')
    print(conf)
    
    refs_all = create_ref_scenes(
        s1_preproc, 
        conf['mode'], conf['pol'],
        conf['ref_start_year'], conf['ref_end_year'],
        ref_doy_start=conf['ref_start_doy'], ref_doy_end=conf['ref_end_doy'], 
        blend_borders=True, blend_borders_nibble_px=conf['export_resolution_m']/25
        )




## Set up a join between summer acquisitions and references, then correct.
# Join criterion for orbits [this implicitly deals with ASC/DESC]
filter_orbits = ee.Filter.equals(leftField='relativeOrbitNumber_start', rightField='orbit') 

# Join criterion for mode
filter_mode = ee.Filter.equals(leftField='instrumentMode', rightField='mode')

# Join criterion for years
if conf['annual_refs'] is True:
    filter_years1 = ee.Filter.equals(leftField='year', rightField='year_start')
    filter_years2 = ee.Filter.equals(leftField='year', rightField='year_end')       
    # Combine the join criteria
    filter_combined = ee.Filter.And(filter_orbits, filter_years1, filter_years2, filter_mode)
else:
    refs = ( refs_all
        .filter(ee.Filter.eq('year_start', conf['ref_start_year']))
        .filter(ee.Filter.eq('year_end',   conf['ref_end_year']))
            )
    filter_combined = ee.Filter.And(filter_orbits, filter_mode)

# Only retain reference scenes which were generated from sufficient scenes
# Threshold=50 is based on Buth et al. (in prep.)
refs_valid = refs.filter(ee.Filter.gt('nscenes', 50))

# Find all summer/autumn/start-of-winter images
all_summer = s1_preproc.filter(ee.Filter.dayOfYear(120, 349))  # 152 to 304 = 1 June to 31 Oct

# Give them a year attribute
def add_year(im):
    year = ee.Date(im.get('system:time_start')).get('year')
    return im.set('year',year)
all_summer = all_summer.map(add_year)

# Apply custom despeckling ...
# This has to happen BEFORE debiasing
# Lee filters require linear scale (_FLOAT), they do not work with log data
fdespeckle = prepare_despeckle(conf['dfilter'], conf['filterpx'], conf['mode'], conf['pol'])
all_summer_despeck = all_summer.map(fdespeckle)
# Convert FLOAT to LOG for the rest of the processing
if 'FLOAT' in conf['collection']:
    print('Converting float (linear) to dB (log)')
    def toDbBands(im):
        return toDb(im, bands=[conf['pol']])
    all_summer = all_summer_despeck.map(toDbBands)
else:
    print('Warning: despeckling took place on dB (log) data. Was this intentional?')
    all_summer = all_summer_despeck
# def toDbBands(im):
#     return toDb(im, bands=[conf['pol']])
# all_summer = all_summer.map(toDbBands)

if debias_method == 'join':
    # Use the Join criteria to match up summer images with a reference scene
    valid_summer_join = ee.Join.saveFirst('reference')
    # This is a FeatureCollection (not an ImageCollection)
    valid_summer = valid_summer_join.apply(all_summer, refs_valid, filter_combined)    
    def debias_pol_band(im):
        return debias_collection_via_join(im, conf['pol'])    
elif debias_method == 'bylist':
    # Make a list of available relative orbit references
    if conf['ref_type'] == 'asset':
        info = refs_valid.getInfo() 
        orbits = [o['properties']['orbit'] for o in info['features']]
        # retain only summer scenes which have a ref scene
        valid_summer = all_summer.filter(ee.Filter.inList('relativeOrbitNumber_start', orbits))
    else: # 'dynamic' case
        valid_summer = all_summer
    def debias_pol_band(im):
        return debias_collection_via_lookup(im, conf['pol'], refs_valid)

s1_debiased = ee.ImageCollection(valid_summer.map(debias_pol_band))

def melt_threshold(im, threshold=-2.1): #was using -3
    # Based on Liang et al. (2021) - however their percentile analysis yields -2.66dB for their study.
    mod = im.select(conf['pol']).lt(threshold) 
    return im.updateMask(mod)

s1_melt = s1_debiased.map(melt_threshold)





