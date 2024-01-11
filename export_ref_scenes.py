"""
Export reference scenes of each orbit as GEE Assets.

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
import argparse
import tomli

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
    

    p = argparse.ArgumentParser(description='Export Sentinel-1 reference scenes as Assets')
    p.add_argument('-restart', type=int, default=None, help='int, number of collection to restart at')
    p.add_argument('config_file', type=str, help='str, path to job configuration TOML file')
    p.add_argument('-debug', dest='debug', default=False, action='store_true', help='Run in debug mode')
    args = p.parse_args()

    with open(args.config_file, "rb") as f:
        conf = tomli.load(f)

    s1 = import_collection(conf['collection'], conf['mode'], conf['pol'], 
        do_remove_dark_borders=False,
        do_mask_floor=False)
    s1 = s1.filterBounds(bounds_gris)
    

    if args.debug:
        s1 = s1.filter(ee.Filter.eq('relativeOrbitNumber_start', 90))
        bounds_test = ee.Geometry.Polygon(
            [[-50.405273, 72.40545],
            [-50.405273, 73.48140],
            [-48.449707, 73.48140],
            [-48.449707, 72.40545],
            [-50.405273, 72.40545]],
            geodesic=False, proj='epsg:4326'
        )
        #PROCESS_REGION = test_gris
        #fn_append = 'TEST'
        PROCESS_REGION = bounds_gris
        fn_append = ''
    else:
        PROCESS_REGION = bounds_gris
        fn_append = ''

    # Do despeckling as a function of export resolution.

    def compute_at_res(im):
        return im.select(conf['pol']).reproject(im.select(conf['pol']).projection(), scale=25) # 50m is 2x 25 IW resol
    
    if conf['dfilter'] is not False:
        fdespeckle = prepare_despeckle(conf['dfilter'], conf['filterpx'], conf['mode'], conf['pol'])
        s1 = s1.map(fdespeckle)
        # Convert FLOAT to LOG for the rest of the processing
        if 'FLOAT' in conf['collection']:
            print('Converting float (linear) to dB (log)')
            def toDbBands(im):
                return toDb(im, bands=[conf['pol']])
            s1 = s1.map(toDbBands)
        else:
            print('Warning: despeckling took place on dB (log) data. Was this intentional?')
        despeckle_suffix = '_{f}{p}'.format(f=conf['dfilter'], p=conf['filterpx'])
    else:
        despeckle_suffix = ''

    ## Force despeckling at high resolution 
    s1 = s1.map(compute_at_res)

    refs = None    
    if conf['annual_refs']:
        print('Producing annual reference scenes')
        for year in range(conf['ref_start_year'], conf['ref_end_year']):
            refs_year = create_ref_scenes(s1.select(conf['pol']), conf['mode'], conf['pol'], year, year,
                ref_doy_start=conf['ref_start_doy'], ref_doy_end=conf['ref_end_doy'])
            if refs is not None:
                refs = refs.merge(refs_year)    
            else:
                refs = refs_year
    else:
        refs = create_ref_scenes(
            s1, conf['mode'], conf['pol'],
            conf['ref_start_year'], conf['ref_end_year'],
            ref_doy_start=conf['ref_start_doy'], ref_doy_end=conf['ref_end_doy'],
            blend_borders=False
            )  

    def export_ref(im, region):
        kws = dict(crs=conf['crs'], scale=conf['export_resolution_m'], maxPixels=1e13)
        fn_prefix = 'ref_{mode}_{pol}_{ron}_{ys}_{ye}_{ds}_{de}_{res}m{d}{test}'
        im_info = im.getInfo()
        fn_prefix = fn_prefix.format(
            mode=im_info['properties']['mode'],
            pol=im_info['bands'][0]['id'],
            ron=im_info['properties']['orbit'],
            ys=im_info['properties']['year_start'], #dt.datetime.fromtimestamp(im_info['properties']['system:time_start']/1000).year, #conf['ref_start_year'],
            ye=im_info['properties']['year_end'], #int(ee.Date(im.get('system:time_end')).get('year').getInfo()), #conf['ref_end_year'],
            ds=conf['ref_start_doy'],
            de=conf['ref_end_doy'],
            res=conf['export_resolution_m'],
            d=despeckle_suffix,
            test=fn_append)
        task = ee.batch.Export.image.toAsset(
            image=im.unmask(-9999),
            region=region,
            description=fn_prefix,
            assetId='{a}/{p}_WITHOUTBORDERCLEAN'.format(a=conf['assets_path'], p=fn_prefix), 
            **kws)
        task.start()
        return task

    print('Calculating number of scenes to export...')
    nrefs = refs.size().getInfo()
    if args.restart is not None:
        offset = args.restart
        print('Starting at {o}'.format(o=offset))
    else:
        offset = 0
        print('{n} reference scenes to export ...'.format(n=nrefs))
    while offset < nrefs:
        print('Batching: {n}-{z}'.format(n=offset, z=(offset + RETR_MAX)))

        img_lst = refs.toList(RETR_MAX, offset=offset)
        imgs_in_lst = img_lst.size().getInfo()
        task_store = []
        count = 0

        for i in range(0, imgs_in_lst):
            image = ee.Image(img_lst.get(i))

            t = export_ref(image.select(conf['pol']), PROCESS_REGION)

            status = t.status()
            print('{n}/{t} {d} : {s}'.format(n=count+1, t=imgs_in_lst, 
                    d=status['description'], s=status['state']))

            task_store.append(t)
            count += 1

        offset += (RETR_MAX + 1)
