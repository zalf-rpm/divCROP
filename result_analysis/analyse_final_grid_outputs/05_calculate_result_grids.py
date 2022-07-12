import os
import glob
import numpy as np
from numpy.core.numeric import outer
from osgeo import gdal
import joblib

WD = "/beegfs/jaenicke/klimertrag_temp/raster/sim-yields/11_final_results"
OUTPUT_FOLDER = f"{WD}/results"

OUTPUT_VARS = ["Yield"]  #"Yield"

SETUPS = {
    # "barleywinterbarley":{
    #     "crop" : "barleywinterbarley",
    #     "rcp26" : [2,3,8,12,13],
    #     "rcp45" : [4,5,7,10,14,15],
    #     "rcp85" : [1,6,9,11,16,17],
    #     "hist_per" : {"start": 1971, "end": 2000},
    #     "futur_per" : {"start": 2031, "end": 2060},
    #     "futur_per2" : {"start": 2070, "end": 2099},
    # },
    # "barleyspringbarley":{
    #     "crop" : "barleyspringbarley",
    #     "rcp26" : [20,21,26,30,31],
    #     "rcp45" : [22,23,25,28,32,33],
    #     "rcp85" : [19,24,27,29,34,35],
    #     "hist_per" : {"start": 1971, "end": 2000},
    #     "futur_per" : {"start": 2031, "end": 2060},
    #     "futur_per2" : {"start": 2070, "end": 2098},
    # },
    # "wheatwinterwheat":{
    #     "crop" : "wheatwinterwheat",
    #     "rcp26" : [38,39,44,48,49],
    #     "rcp45" : [40,41,43,46,50,51],
    #     "rcp85" : [37,42,45,47,52,53],
    #     "hist_per" : {"start": 1971, "end": 2000},
    #     "futur_per" : {"start": 2031, "end": 2060},
    #     "futur_per2" : {"start": 2070, "end": 2099},
    # },
    # "ryewinterrye":{
    #     "crop" : "ryewinterrye",
    #     "rcp26" : [56,57,62,66,67],
    #     "rcp45" : [58,59,61,64,68,69],
    #     "rcp85" : [55,60,63,65,70,71],
    #     "hist_per" : {"start": 1971, "end": 2000},
    #     "futur_per" : {"start": 2031, "end": 2060},
    #     "futur_per2" : {"start": 2070, "end": 2099},
    # },
    # "rapewinterrape":{
    #     "crop" : "rapewinterrape",
    #     "rcp26" : [74,75,80,84,85],
    #     "rcp45" : [76,77,79,82,86,87],
    #     "rcp85" : [73,78,81,83,88,89],
    #     "hist_per" : {"start": 1971, "end": 2000},
    #     "futur_per" : {"start": 2031, "end": 2060},
    #     "futur_per2" : {"start": 2070, "end": 2099},
    # },
    "maizesilagemaize":{
        "crop" : "maizesilagemaize",
        "rcp26" : [92,93,98,102,103],
        "rcp45" : [94,95,97,100,104,105],
        "rcp85" : [91,96,99,101,106,107],
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2098},
    },
    "sugarbeet":{
        "crop" : "sugarbeet",
        "rcp26" : [110,111,116,120,121],
        "rcp45" : [112,113,115,118,122,123],
        "rcp85" : [109,114,117,119,124,125],
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2098},
    },
    "potatomoderatelyearlypotato":{
        "crop" : "potatomoderatelyearlypotato",
        "rcp26" : [128,129,134,138,139],
        "rcp45" : [130,131,133,136,140,141],
        "rcp85" : [127,132,135,137,142,143],
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2098},
    },
    }

def create_folder(directory):
    """
    Tries to create a folder at the specified location. Path should already exist (excluding the new folder). 
    If folder already exists, nothing will happen.
    :param directory: Path including new folder.
    :return: Creates a new folder at the specified location.
    """

    import os
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' + directory )

def write_array_to_raster(in_array,
                          out_path,
                          gt,
                          pr,
                          no_data_value,
                          type_code=None,
                          options=['COMPRESS=DEFLATE', 'PREDICTOR=1'],
                          driver='GTiff'):
    """
    Writes an array to a tiff-raster. If no type code of output is given, it will be extracted from the input array.
    As default a deflate compression is used, but can be specified by the user.
    :param in_array: Input array
    :param out_path: Path of output raster
    :param gt: GeoTransfrom of output raster
    :param pr: Projection of output raster
    :param no_data_value: Value that should be recognized as the no data value
    :return: Writes an array to a raster file on the disc.
    """

    from osgeo import gdal
    from osgeo import gdal_array

    if type_code == None:
        type_code = gdal_array.NumericTypeCodeToGDALTypeCode(in_array.dtype)

    if len(in_array.shape) == 3:
        nbands_out = in_array.shape[0]
        x_res = in_array.shape[2]
        y_res = in_array.shape[1]

        out_ras = gdal.GetDriverByName(driver).Create(out_path, x_res, y_res, nbands_out, type_code, options=options)
        out_ras.SetGeoTransform(gt)
        out_ras.SetProjection(pr)

        for b in range(0, nbands_out):
            band = out_ras.GetRasterBand(b + 1)
            arr_out = in_array[b, :, :]
            band.WriteArray(arr_out)
            band.SetNoDataValue(no_data_value)
            band.FlushCache()

        del (out_ras)

    if len(in_array.shape) == 2:
        nbands_out = 1
        x_res = in_array.shape[1]
        y_res = in_array.shape[0]

        out_ras = gdal.GetDriverByName(driver).Create(out_path, x_res, y_res, nbands_out, type_code, options=options)
        out_ras.SetGeoTransform(gt)
        out_ras.SetProjection(pr)

        band = out_ras.GetRasterBand(1)
        band.WriteArray(in_array)
        band.SetNoDataValue(no_data_value)
        band.FlushCache()

        del (out_ras)

def calculate_statistics(file_lst, statistic, out_pth):

    if statistic not in ['mean', 'std', 'cov']:
        print("Please choose between 'mean', 'std' and 'cov' for calculating the \
            mean, the standard deviation and the coefficient of variation, respectively.")
        return

    ## check if all files exist
    for file in file_lst:
        if not os.path.exists(file):
            print(file, "from file list does not exist!")
            return

    ## load data
    ras_lst = [gdal.Open(file) for file in file_lst]
    arr_lst = [ras.ReadAsArray() for ras in ras_lst]
    ndv_lst = [ras.GetRasterBand(1).GetNoDataValue() for ras in ras_lst]

    ## mask arrays
    arr_m_lst = []
    for i, arr in enumerate(arr_lst):
        ndv = ndv_lst[i]

        ## Create a mask that sets all no data values to 0
        ndv_mask = np.where(arr == ndv, 0, 1)

        ## Mask yield array with no data mask
        arr_m = np.ma.masked_where(ndv_mask == 0, arr)
        arr_m_lst.append(arr_m.copy())

    if statistic == 'mean':
        ## Average yields over all years
        res_arr = np.ma.mean(arr_m_lst, axis=0)
    elif statistic == 'std':
        ## Calculate standard deviation
        res_arr = np.ma.std(arr_m_lst, axis=0)
    elif statistic == 'cov':
        ## Calculate coefficient of variation
        std_arr = np.ma.std(arr_m_lst, axis=0)
        avg_arr = np.ma.mean(arr_m_lst, axis=0)
        res_arr = std_arr / avg_arr

    gt = ras_lst[0].GetGeoTransform()
    pr = ras_lst[0].GetProjection()

    res_arr = res_arr.filled(ndv)

    write_array_to_raster(res_arr,
                        out_path=out_pth,
                        gt=gt,
                        pr=pr,
                        no_data_value=ndv,
                        type_code=None,
                        options=['COMPRESS=DEFLATE', 'PREDICTOR=1'],
                        driver='GTiff')

    print(f"Calculated {statistic} of file list. File: {out_pth}.")

def calculate_difference_between_rasters(file_minuend, file_subtrahend, no_data_value, out_pth, percentage_change=False):
    """
    Calculates difference between 1-band arrays (arr_minuend - arr_subtrahend). Output can be absolute change or %-change.
    """
    ras_minuend = gdal.Open(file_minuend)
    ras_subtrahend = gdal.Open(file_subtrahend)

    arr_minuend = ras_minuend.ReadAsArray()
    arr_subtrahend = ras_subtrahend.ReadAsArray()

    ndv_minu = ras_minuend.GetRasterBand(1).GetNoDataValue()
    ndv_subt = ras_subtrahend.GetRasterBand(1).GetNoDataValue()

    if arr_subtrahend.ndim > 2:
        print("Array 1 has more than 2 dimensions! Only two dimensions allowed!")
        return

    if arr_minuend.ndim > 2:
        print("Array 2 has more than 2 dimensions!  Only two dimensions allowed!")
        return

    if arr_subtrahend.shape != arr_minuend.shape:
        print("Provided array have different shapes. Please provide array with same shapes and only 1 band.")
        return

    arr_subtrahend_m = np.where(arr_subtrahend == ndv_subt, 0, 1)
    arr_minuend_m = np.where(arr_minuend == ndv_minu, 0, 1)

    mask = arr_subtrahend_m * arr_minuend_m

    diff_arr = arr_minuend - arr_subtrahend

    if percentage_change:
        diff_arr = np.divide(diff_arr, arr_subtrahend) * 100

    diff_arr = np.where(mask == 1, diff_arr, no_data_value)

    gt = ras_minuend.GetGeoTransform()
    pr = ras_minuend.GetProjection()

    write_array_to_raster(diff_arr,
                        out_path=out_pth,
                        gt=gt,
                        pr=pr,
                        no_data_value=no_data_value,
                        type_code=None,
                        options=['COMPRESS=DEFLATE', 'PREDICTOR=1'],
                        driver='GTiff')

    return diff_arr

def calculate_yearly_statistics(wd, setup, rcp, period, output_var, out_folder, statistic):

    create_folder(out_folder)

    run_ids = setup[rcp]
    crop = setup["crop"]

    start = setup[period]["start"]
    end = setup[period]["end"]

    for year in range(start, end+1):
        file_lst = []
        for run_id in run_ids:
            search_term = f"{wd}/{run_id}/{crop}_{output_var}_{year}*.asc"
            file_pth = glob.glob(search_term)
            if len(file_pth) == 1:
                file_pth = file_pth[0]
                file_lst.append(file_pth)
            elif len(file_pth) > 1:
                print(f"Two files found for RunID {run_id} and variable '{output_var}'. There should only be one!")
            else:
                print(f"No file found for RunID {run_id} and variable '{output_var}'.")
        
        if len(file_lst) < 1:
            print(f"No files found. {crop}, {output_var}, {rcp}, {year}")
            continue
        else:
            print(f"For {crop}, {output_var}, {rcp}, {year} there are {len(file_lst)} files.")
        
        out_pth = f"{out_folder}/{crop}_{rcp}_{output_var}_{year}_{statistic}.tiff"
        calculate_statistics(file_lst, statistic, out_pth)

def calculate_period_statistics(wd, setup, rcp, period, output_var, out_folder, statistic):

    create_folder(out_folder)
    crop = setup["crop"]

    start = setup[period]["start"]
    end = setup[period]["end"]

    file_lst = []
    for year in range(start, end+1):
        file_pth = f"{wd}/{crop}_{rcp}_{output_var}_{year}_mean.tiff"
        file_lst.append(file_pth)
        
    if len(file_lst) != (end-start):
        print(f"There are {len(file_lst)} for {end-start} years.")
    
    out_pth = f"{out_folder}/{crop}_{rcp}_{output_var}_{start}-{end}_{statistic}.tiff"
    calculate_statistics(file_lst, statistic, out_pth)

def calculate_change(wd, setup, rcp, period1, period2, output_var, out_folder):

    create_folder(out_folder)

    crop = setup["crop"]

    start_fut = setup[period2]["start"]
    end_fut = setup[period2]["end"]
    start_hist = setup[period1]["start"]
    end_hist = setup[period1]["end"]

    file_future = f"{wd}/{crop}_{rcp}_{output_var}_{start_fut}-{end_fut}_mean.tiff"
    file_hist = f"{wd}/{crop}_{rcp}_{output_var}_{start_hist}-{end_hist}_mean.tiff"
    out_pth = f"{out_folder}/{crop}_{rcp}_{output_var}_diff_{start_fut}-{end_fut}_{start_hist}-{end_hist}.tiff"    
    calculate_difference_between_rasters(file_minuend=file_future, file_subtrahend=file_hist, no_data_value=0, out_pth=out_pth, percentage_change=False)

    out_pth = f"{out_folder}/{crop}_{rcp}_{output_var}_diff_{start_fut}-{end_fut}_{start_hist}-{end_hist}_perc.tiff"    
    calculate_difference_between_rasters(file_minuend=file_future, file_subtrahend=file_hist, no_data_value=0, out_pth=out_pth, percentage_change=True)

def work_func(sname):
    setup = SETUPS[sname]
    print(setup)
    for output_var in OUTPUT_VARS:
        ## We have to use a mask

        ## 1. Calculate the average per year across the different rcms
       
        out_yearly = f"{OUTPUT_FOLDER}/yearly_averages/{sname}"
        for rcp in ["rcp26", "rcp45", "rcp85"]: #
            for period in ["hist_per", "futur_per", "futur_per2"]:
                calculate_yearly_statistics(wd=WD, setup=setup, rcp=rcp, period=period, output_var=output_var, out_folder=out_yearly, statistic="mean")

        ## 2. Calculate the average, std and the coeff of var. per rcp across the time periods
        out_periods = f"{OUTPUT_FOLDER}/period_statistics/{sname}"
        for rcp in ["rcp26", "rcp45", "rcp85"]:
            for period in [ "futur_per2"]:  #"hist_per", "futur_per",
                calculate_period_statistics(wd=out_yearly, setup=setup, rcp=rcp, period=period, output_var=output_var, out_folder=out_periods, statistic="mean")
                calculate_period_statistics(wd=out_yearly, setup=setup, rcp=rcp, period=period, output_var=output_var, out_folder=out_periods, statistic="std")
                calculate_period_statistics(wd=out_yearly, setup=setup, rcp=rcp, period=period, output_var=output_var, out_folder=out_periods, statistic="cov")
            
        ## 3. Calculate the difference between historical averages and future averages
        out_change= f"{OUTPUT_FOLDER}/change_maps/{sname}"
        for rcp in ["rcp26", "rcp45", "rcp85"]:
            calculate_change(wd=out_periods, setup=setup, rcp=rcp, period1="hist_per", period2="futur_per", output_var=output_var, out_folder=out_change)
            calculate_change(wd=out_periods, setup=setup, rcp=rcp, period1="hist_per", period2="futur_per2", output_var=output_var, out_folder=out_change)
    

def main():
    for output_var in OUTPUT_VARS:
        for sname in SETUPS:

            ## We have to use a mask

            ## 1. Calculate the average per year across the different rcms
            setup = SETUPS[sname]
            out_yearly = f"{OUTPUT_FOLDER}/yearly_averages/{sname}"
            for rcp in ["rcp26", "rcp45", "rcp85"]: #
                for period in ["hist_per", "futur_per", "futur_per2"]:
                    calculate_yearly_statistics(wd=WD, setup=setup, rcp=rcp, period=period, output_var=output_var, out_folder=out_yearly, statistic="mean")

            ## 2. Calculate the average, std and the coeff of var. per rcp across the time periods
            out_periods = f"{OUTPUT_FOLDER}/period_statistics/{sname}"
            for rcp in ["rcp26", "rcp45", "rcp85"]:
                for period in ["hist_per", "futur_per", "futur_per2"]:
                    calculate_period_statistics(wd=out_yearly, setup=setup, rcp=rcp, period=period, output_var=output_var, out_folder=out_periods, statistic="mean")
                    calculate_period_statistics(wd=out_yearly, setup=setup, rcp=rcp, period=period, output_var=output_var, out_folder=out_periods, statistic="std")
                    calculate_period_statistics(wd=out_yearly, setup=setup, rcp=rcp, period=period, output_var=output_var, out_folder=out_periods, statistic="cov")
                
            ## 3. Calculate the difference between historical averages and future averages
            out_change= f"{OUTPUT_FOLDER}/change_maps/{sname}"
            for rcp in ["rcp26", "rcp45", "rcp85"]:
                calculate_change(wd=out_periods, setup=setup, rcp=rcp, period1="hist_per", period2="futur_per", output_var=output_var, out_folder=out_change)
                calculate_change(wd=out_periods, setup=setup, rcp=rcp, period1="hist_per", period2="futur_per2", output_var=output_var, out_folder=out_change)
        
if __name__ == '__main__':
    main()
    # joblib.Parallel(n_jobs=9)(joblib.delayed(work_func)(sname) for sname in SETUPS)
