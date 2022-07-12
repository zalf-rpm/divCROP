import matplotlib.pyplot as plt
import pandas as pd
import glob
from functools import reduce
import os
import numpy as np
from sklearn import metrics
import math
from osgeo import gdal
import seaborn as sns
from matplotlib.lines import Line2D

RUN_IDS = [18, 36, 54, 72, 90, 108, 126, 144]  # 

## Working directory
WD = "/beegfs/jaenicke/klimertrag/result_analysis"

## Reference values for phenology and yields
OBS_PTH = WD + "/data/reference_tables/phenology/PH_Jahresmelder_Landwirtschaft_{0}_hist_filt_years.csv"  # DWD Crop name
OBS_PTH_GER = f"{WD}/data/reference_tables/yields/Ertraege_Deutschland_10_Frucharten_1999-2020_detrended.csv"
YIELD_OBS_COL = 'yield_obs_detr'  # column in yield reference data

## Mask to be used for masking the results. Leave empty string when no mask should be applied.
MASK_PTH = "/beegfs/jaenicke/klimertrag/result_analysis/data/crop_masks/CM_2017-2019_{0}_1000m_25832_q3.asc"  # Crop abbreveation
MASK_PTH = "/beegfs/jaenicke/klimertrag/result_analysis/data/crop_masks/CTM_17-19_mask_1000m_25832.asc"

## Folder with outputs
GRID_FOLDER = "/beegfs/jaenicke/klimertrag_temp/raster/sim-yields/11_final_results/{0}"  # RUN_ID

## Output folder for figures
FIGURE_FOLDER = "/beegfs/jaenicke/klimertrag_temp/figures/11_final_results/phenology_validation/{0}"  # RUN_ID

## Events as provided in output csvs
EVENTS = ["Sowing", "cereal-stem-elongation", "Stage-2", "Stage-3", "Stage-4","Stage-5", "Stage-6", "Stage-7", "Harvest"]


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)

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

def d_modified(cal, obs, order=1):
    import warnings
    """
    Modify Index of Agreement (Willmott et al., 1984) range from 0.0 to 1.0
    and the closer to 1 the better the performance of the model.
    :param cal: (N,) array_like of calculated values
    :param obs: (N,) array_like of observed values
    :param order: exponent to be used in the computation. Default is 1
    :return: Modified Index of Agreement between 'cal' and 'obs'
    .. math::
        r = 1 - \\frac{\\sum (O_i - C_i)^n}
                      {\\sum (\\vert C_i - m_O \\vert + \\vert O_i - m_O \\vert)^n}
    where :math:`m_O` is the mean of the vector :math:`O` of observed values.

    Found @ https://github.com/SimonDelmas/goodness_of_fit/blob/master/goodness_of_fit/__init__.py
    """

    obs_mean = np.nanmean(obs)
    denominator = np.nansum(
        np.power(np.abs(cal - obs_mean) + np.abs(obs - obs_mean), order)
    )
    if denominator == 0.0:
        warnings.warn("Index of agreement potential error is null! Return NaN.")
        return np.nan

    nominator = np.nansum(np.abs(np.power(obs - cal, order)))
    return 1 - (nominator / denominator)

def monica_grids_results_to_df(folder, event, mask_pth=None):

    event_to_var_name = {
        "Sowing": "sdoy",
        "Stage-2": "s2doy",
        "cereal-stem-elongation": "sedoy",
        "Stage-3": "s3doy",
        "Stage-4": "s4doy",
        "Stage-5": "s5doy",
        "Stage-6": "s6doy",
        "Stage-7": "s7doy",
        "Harvest": "hdoy"
    }

    var = event_to_var_name[event]
    pth_lst = glob.glob(f"{folder}/*{var}*.asc")

    year_lst = []
    arr_lst = []

    for pth in pth_lst:
        bname = os.path.basename(pth)
        year = bname.split('_')[2]
        year_lst.append(year)

        ras = gdal.Open(pth)
        arr = ras.ReadAsArray()
        if mask_pth:
            ras_mask = gdal.Open(mask_pth)
            mask = ras_mask.ReadAsArray()
            arr = np.where(mask == 1, arr, -9999)
        arr = arr.flatten()
        arr = arr[arr != -9999]
        arr_lst.append(arr)
    
    if len(arr_lst) == 0:
        df = pd.DataFrame(columns=["empty"])
        return df

    min_x = min([arr.shape[0] for arr in arr_lst])
    arr_lst = [arr[:min_x] for arr in arr_lst]
    final_arr = np.vstack(arr_lst)
    final_arr = final_arr.T

    df = pd.DataFrame(final_arr)
    df.columns = year_lst

    return df

def plot_simulations_vs_observations(df_plt, out_pth):

    plt.figure(figsize=(20, 10))
    ax = sns.boxplot(y='value', x='Year',
                hue='type',
                data=df_plt,
                palette="colorblind",
                fliersize=0.1)
    # ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
    plt.grid('minor')
    plt.savefig(out_pth)
    plt.close()
    print('Plotting done!', out_pth)

def main():

    for run_id in RUN_IDS:

        grid_folder = GRID_FOLDER.format(run_id)

        crop_name = glob.glob(f"{grid_folder}/*.asc")
        crop_name = [os.path.basename(item).split('_')[0] for item in crop_name]
        crop_name = list(set(crop_name))[0]

        english_crop_name_to_dwd_crop_name = {
            "barleywinterbarley": "Wintergerste",
            "barleyspringbarley": "Sommergerste",
            "wheatwinterwheat": "Winterweizen",
            "ryewinterrye": "Winterroggen",
            "rapewinterrape": "Winterraps",
            "maizesilagemaize": "Mais",
            "sugarbeet": "Zucker-Ruebe"
        }

        if not crop_name in english_crop_name_to_dwd_crop_name:
            print(f"{crop_name} is not in DWD crop names. Skipping.")
            continue

        dwd_crop_name = english_crop_name_to_dwd_crop_name[crop_name]

        for event in EVENTS:
            df_obs = pd.read_csv(OBS_PTH.format(dwd_crop_name))

            monica_to_dwd_stages = {
                "Wintergerste": {
                    "Sowing": 10,
                    "Stage-2": 12,
                    "cereal-stem-elongation": 15,
                    "Stage-3": 15,
                    "Stage-4": 18,
                    "Stage-6": 21,
                    "Harvest": 24
                },
                "Sommergerste": {
                    "Sowing": 10,
                    "Stage-2": 12,
                    "cereal-stem-elongation": 15,
                    "Stage-3": 15,
                    "Stage-4": 18,
                    "Stage-6": 21,
                    "Harvest": 24
                },
                "Winterweizen": {
                    "Sowing": 10,
                    "Stage-2": 12,
                    "cereal-stem-elongation": 15,
                    "Stage-3": 15,
                    "Stage-4": 18,
                    "Stage-5": 19,
                    "Stage-6": 21,
                    "Harvest": 24
                },
                "Winterroggen": {
                    "Sowing": 10,
                    "Stage-2": 12,
                    "cereal-stem-elongation": 15,
                    "Stage-3": 15,
                    "Stage-4": 18,
                    "Stage-5": 5,
                    "Stage-6": 21,  #6,
                    "Harvest": 24
                },
                "Winterraps": {
                    "Sowing": 10,
                    "Stage-2": 12,
                    "cereal-stem-elongation": 67, ## Das stimmt sehr wahrschenlich nicht
                    "Stage-3": 67, ##14,  ## oder 67
                    "Stage-4": 17,
                    "Stage-5": 5,
                    "Stage-6": 22,  #6,  #
                    "Harvest": 24
                },
                "Mais": {
                    "Sowing": 10,
                    "Stage-2": 12,
                    "cereal-stem-elongation": 67,  ## Das stimmt sehr wahrschenlich nicht
                    "Stage-3": 67,  ##15 oder 67
                    "Stage-4": 65,
                    "Stage-5": 5,
                    "Stage-6": 19,  #
                    "Stage-7": 20,  ## oder 21
                    "Harvest": 24
                },
                "Zucker-Ruebe": {
                    "Sowing": 10,
                    "Stage-2": 12,
                    "Stage-3": 13,
                    "Stage-4": 13,
                    "Harvest": 24
                }
            }

            if not event in monica_to_dwd_stages[dwd_crop_name]:
                print(f"{event} is not in DWD phases for {dwd_crop_name}. Skipping.")
                continue

            phase_id = monica_to_dwd_stages[dwd_crop_name][event]
            df_obs = df_obs.loc[df_obs['Phase_id'] == phase_id].copy()

            df_obs = pd.melt(df_obs[['Referenzjahr', 'Jultag']], id_vars='Referenzjahr', value_vars='Jultag',
                                    var_name='type', value_name='value')
            df_obs.columns = ['Year', 'type', 'value']
            df_obs['type'] = 'obs'

            ### Test against grid simulations 
            if os.path.exists(grid_folder):
                crop_name_to_abbrevation = {
                    "Wintergerste" : "WB",
                    "Sommergerste" : "SB",
                    "Winterweizen" : "WW",
                    "Mais": "SM",
                    "Zucker-Ruebe" : "SU",
                    "Winterroggen": "WRye",
                    "Kartoffel": "PO",
                    "Winterraps": "WR"
                }

                if not dwd_crop_name in crop_name_to_abbrevation:
                    print(f"{dwd_crop_name} is not in Crop abbrevations. Skipping.")
                    continue

                crop_abbr = crop_name_to_abbrevation[dwd_crop_name]
                df_sim = monica_grids_results_to_df(grid_folder, event, mask_pth=MASK_PTH.format(crop_abbr))
                if df_sim.empty():
                    print(f"No grids found for {event} of {dwd_crop_name}")
                    continue
                
                df_sim = pd.melt(df_sim, value_vars=df_sim.columns.tolist(), var_name='Year', value_name='value')
                df_sim['type'] = 'sim_grid'
                df_sim = df_sim[['Year', 'type', 'value']]

                df_plt = pd.concat([df_sim, df_obs])
                df_plt['Year'] = df_plt['Year'].astype(int)
                df_plt = df_plt.loc[(df_plt['Year'] > 1998) & (df_plt['Year'] < 2020)].copy()

                out_pth = f"{FIGURE_FOLDER.format(run_id)}/{dwd_crop_name}_{event}_grid_simulations_boxplot_dwd_pheno.png"
                create_folder(os.path.dirname(out_pth))
                plot_simulations_vs_observations(df_plt, out_pth)
   

if __name__ == '__main__':
    main()