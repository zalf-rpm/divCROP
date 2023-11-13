#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import os
import json
import numpy as np
from pyproj import CRS, Transformer
import sqlite3
import sys
import time as times

import soil_io3
import monica_run_lib as Mrunlib

PATHS = {
    "hpc-local": {
        "path-to-climate-dir": "/beegfs/common/data/climate/",
        "monica-path-to-climate-dir": "/beegfs/common/data/climate/",
        "path-to-data-dir": "./data/",  
        "out-folder": "./out/",
    },
}

DATA_SOIL_DB = "germany/buek200.sqlite"
DATA_GRID_LAND_USE = "germany/landuse_1000_31469_gk5.asc"
DATA_GRID_SOIL = "germany/buek200_1000_25832_etrs89-utm32n.asc"
DATA_GRID_CROPS = "germany/crops-all2017-2019_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/germany-complete_1000_25832_etrs89-utm32n.asc"
TEMPLATE_PATH_LATLON = "{path_to_climate_dir}/latlon-to-rowcol.json"
TEMPLATE_PATH_CLIMATE_CSV = "{gcm}/{rcm}/{scenario}/{ensmem}/{version}/row-{crow}/col-{ccol}.csv"


def run_climate_analyzer():
    """main"""

    config = {
        "mode": "hpc-local",  
        "path_to_dem_grid": "",
        "setups-file": "setups_climate.csv",
        "merge-setups": "[1]",
    }

    # read commandline args only if script is invoked directly from commandline
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    print("config:", config)
    print("current folder", os.getcwd())

    # select paths 
    paths = PATHS[config["mode"]]
    # open soil db connection
    soil_db_con = sqlite3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB)

    # read setup from csv file
    setups = Mrunlib.read_sim_setups(config["setups-file"])
    run_setups = json.loads(config["merge-setups"])
    print("read sim setups: ", config["setups-file"])

    # transforms geospatial coordinates from one coordinate reference system to another
    # transform wgs84 into gk5
    soil_crs_to_x_transformers = {}
    wgs84_crs = CRS.from_epsg(4326)

    # Load grids
    ## note numpy is able to load from a compressed file, ending with .gz or .bz2

    # soil data
    path_to_soil_grid = paths["path-to-data-dir"] + DATA_GRID_SOIL
    print(path_to_soil_grid)
    soil_epsg_code = int(path_to_soil_grid.split("/")[-1].split("_")[2])
    soil_crs = CRS.from_epsg(soil_epsg_code)
    if wgs84_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[wgs84_crs] = Transformer.from_crs(soil_crs, wgs84_crs, always_xy=True)
    soil_metadata, header = Mrunlib.read_header(path_to_soil_grid)
    soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    print("read: ", path_to_soil_grid)

    scols = int(soil_metadata["ncols"])
    srows = int(soil_metadata["nrows"])
    scellsize = int(soil_metadata["cellsize"])
    xllcorner = int(soil_metadata["xllcorner"])
    yllcorner = int(soil_metadata["yllcorner"])
    nodata_value = int(soil_metadata["nodata_value"])

    # land use data
    path_to_landuse_grid = paths["path-to-data-dir"] + DATA_GRID_LAND_USE
    landuse_epsg_code = int(path_to_landuse_grid.split("/")[-1].split("_")[2])
    landuse_crs = CRS.from_epsg(landuse_epsg_code)
    if landuse_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[landuse_crs] = Transformer.from_crs(soil_crs, landuse_crs, always_xy=True)
    landuse_meta, _ = Mrunlib.read_header(path_to_landuse_grid)
    landuse_grid = np.loadtxt(path_to_landuse_grid, dtype=int, skiprows=6)
    landuse_interpolate = Mrunlib.create_ascii_grid_interpolator(landuse_grid, landuse_meta)
    print("read: ", path_to_landuse_grid)

    # crop mask data
    path_to_crop_grid = paths["path-to-data-dir"] + DATA_GRID_CROPS
    crop_epsg_code = int(path_to_crop_grid.split("/")[-1].split("_")[2])
    crop_crs = CRS.from_epsg(crop_epsg_code)
    if crop_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[crop_crs] = Transformer.from_crs(soil_crs, crop_crs, always_xy=True)
    crop_grid = np.loadtxt(path_to_crop_grid, dtype=int, skiprows=6)
    print("read: ", path_to_crop_grid)

    # create output folder if not exists
    os.makedirs(paths["out-folder"], exist_ok=True)

    # group setups by merge-id
    setups_by_merge_id = dict()
    for _, merge_setup_id in enumerate(run_setups):
        for setup_id in setups:
            mergeID =  int(setups[setup_id]["merge"])
            if mergeID == merge_setup_id:
                if mergeID not in setups_by_merge_id:
                    setups_by_merge_id[mergeID] = []
                setups_by_merge_id[mergeID].append(setup_id)

    # run calculations for each setup
    for mergeID in setups_by_merge_id:

        # create grid for storing results
        results_grid = np.zeros((srows, scols), dtype=float)

        for setup_id in setups_by_merge_id[mergeID]:

            setup = setups[setup_id]
            gcm = setup["gcm"]
            rcm = setup["rcm"]
            scenario = setup["scenario"]
            ensmem = setup["ensmem"]
            version = setup["version"]

            start = setup["start_date"]
            start = times.strptime(start, "%Y-%m-%d")
            end = setup["end_date"]
            end = times.strptime(end, "%Y-%m-%d")
            years = end.tm_year - start.tm_year + 1
            temp = setup["temp_threshold"]
            temp = float(temp)            

            use_cmip6_climate_data = (not setup["climate_path_to_latlon_file"]
                                    or len(setup["climate_path_to_latlon_file"]) == 0)

            cdict = climate_data_interpolator = None
            if not use_cmip6_climate_data:
                cdict = {}
                # path to latlon-to-rowcol.json
                path = TEMPLATE_PATH_LATLON.format(path_to_climate_dir=paths["path-to-climate-dir"] +
                                                                    setup["climate_path_to_latlon_file"] + "/")
                climate_data_interpolator = \
                    Mrunlib.create_climate_geoGrid_interpolator_from_json_file(path, wgs84_crs, soil_crs, cdict)
                print("created climate_data to gk5 interpolator: ", path)

            # unknown_soil_ids = set()
            soil_id_cache = {}
            print("All Rows x Cols: " + str(srows) + "x" + str(scols))

            # prevent double calculation of climate data
            previous_climate_path = []
            previous_value = 0
            for srow in range(0, srows):
                print(srow, )

                for scol in range(0, scols):
                    soil_id = int(soil_grid[srow, scol])
                    if soil_id == nodata_value:
                        results_grid[srow, scol] = nodata_value
                        continue

                    # get coordinate of clostest climate element of real soil-cell
                    sh = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                    sr = xllcorner + (scellsize / 2) + scol * scellsize
                    # inter = crow/ccol encoded into integer

                    if use_cmip6_climate_data:
                        c_lon_0 = -179.75
                        c_lat_0 = +89.25
                        c_resolution = 0.5

                        clon, clat = soil_crs_to_x_transformers[wgs84_crs].transform(sr, sh)

                        ccol = int((clon - c_lon_0) / c_resolution)
                        crow = int((c_lat_0 - clat) / c_resolution)

                    else:  # all other climate data will use the interpolator
                        crow, ccol = climate_data_interpolator(sr, sh)
                        clat, _ = cdict[(crow, ccol)]

                    crop_grid_id = int(crop_grid[srow, scol])
                    # mask crop_grid_id set no crop area to no data
                    if crop_grid_id != 1:
                        # set nodata
                        # on position srow, scol
                        results_grid[srow, scol] = nodata_value
                        continue



                    if soil_id in soil_id_cache:
                        soil_profile = soil_id_cache[soil_id]
                    else:
                        soil_profile = soil_io3.soil_parameters(soil_db_con, soil_id)
                        soil_id_cache[soil_id] = soil_profile

                    if len(soil_profile) == 0:
                        # set nodata
                        # on position srow, scol
                        results_grid[srow, scol] = nodata_value
                        continue

                    # check if current grid cell is used for agriculture                
                    if setup["landcover"]: 
                        tcoords = {}
                        if landuse_crs not in tcoords:
                            tcoords[landuse_crs] = soil_crs_to_x_transformers[landuse_crs].transform(sr, sh)
                        lur, luh = tcoords[landuse_crs]
                        landuse_id = landuse_interpolate(lur, luh)
                        if landuse_id not in [2, 3, 4]:
                            # set nodata
                            # on position srow, scol
                            results_grid[srow, scol] = nodata_value
                            continue


                    repl_map = {
                        "gcm": gcm, "rcm": rcm, "scenario": scenario, "ensmem": ensmem, "version": version,
                        "crow": str(crow), "ccol": str(ccol)
                    }
                    for key in list(repl_map.keys()):
                        if key not in setup["climate_path_to_csvs"]:
                            repl_map.pop(key)

                    subpath_to_csv = setup["climate_path_to_csvs"].format_map(repl_map)
                    pathToClimateCSV = [paths["monica-path-to-climate-dir"] + subpath_to_csv]
                    if setup["incl_hist"]:
                        if rcm[:3] == "UHO":
                            hist_subpath_to_csv = setup["climate_path_to_csvs"].format_map(
                                dict(repl_map, rcm="CLMcom-CCLM4-8-17", scenario="historical"))
                            pathToClimateCSV.insert(0, paths["monica-path-to-climate-dir"] + hist_subpath_to_csv)

                        elif rcm[:3] == "SMH":
                            hist_subpath_to_csv = setup["climate_path_to_csvs"].format_map(
                                dict(repl_map, rcm="CLMcom-CCLM4-8-17", scenario="historical"))
                            pathToClimateCSV.insert(0, paths["monica-path-to-climate-dir"] + hist_subpath_to_csv)

                        hist_subpath_to_csv = setup["climate_path_to_csvs"].format_map(
                            dict(repl_map, scenario="historical"))
                        pathToClimateCSV.insert(0, paths["monica-path-to-climate-dir"] + hist_subpath_to_csv)
                    #print("pathToClimateCSV:", pathToClimateCSV)


                    # read climate data from csv file
                    # count if temperature is > 30°C
                    # on position srow, scol
                    # us "crow": int(crow), "ccol": int(ccol),
                    if pathToClimateCSV != previous_climate_path:
                        val = count_days_with_temp_above(pathToClimateCSV, temp, start, end, years)
                        previous_value = val
                    else :
                        val = previous_value
                    results_grid[srow, scol] = results_grid[srow, scol] + val
                    previous_climate_path = pathToClimateCSV

        # devide by number of setups to get mean
        numSetups = len(setups_by_merge_id[mergeID])
        # only divide values that are not nodata
        results_grid[results_grid != nodata_value] = results_grid[results_grid != nodata_value] / numSetups

        # save results to file
        path_to_results_grid = paths["out-folder"] + "results_grid_merge_" + str(mergeID) + ".asc"
        np.savetxt(path_to_results_grid, results_grid, fmt="%1.2f", header=header, comments="")


def count_days_with_temp_above(pathToClimateCSV, temp, start, end, years) :
    #TBD read climate data and calculate average number of days where temperature > threshold 
    #print("pathToClimateCSV:", pathToClimateCSV)
    countedDays = 0
    for path in pathToClimateCSV:
        if os.path.exists(path):
            linesToSkip = 2   
            with open(path, "r") as f:
                lines = f.readlines()
                for line in lines:
                    if linesToSkip > 0:
                        linesToSkip = linesToSkip - 1
                        continue
                
                    token = line.split(",")
                    date = token[0]
                    # convert to date
                    date = times.strptime(date, "%Y-%m-%d")

                    tmax = token[3]
                    if date >= start and date <= end:
                        if float(tmax) > temp:
                            countedDays = countedDays + 1
        else:
            print("file not found:", path)

    # devide by number of years to get mean
    avgDays =  countedDays / years
    return  avgDays


if __name__ == "__main__":
    run_climate_analyzer()
