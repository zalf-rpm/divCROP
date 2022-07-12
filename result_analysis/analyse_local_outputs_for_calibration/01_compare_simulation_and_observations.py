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

DISTRICT_IDS = [1]
RUN_ID = 305
CROP = 'Sommergerste'  #"Winterraps" #"Sommergerste" #"Zucker-Ruebe" #"Mais"

## Working directory
WD = r"D:\projects\KlimErtrag\klimertrag\result_analysis"

## Folder with reference values for phenology and yields
OBS_PTH = fr"{WD}\data\reference_tables\phenology\PH_Jahresmelder_Landwirtschaft_{CROP}_hist_filt_years.csv"
OBS_PTH_GER = rf"{WD}\data\reference_tables\yields\Ertraege_Deutschland_10_Frucharten_1999-2020_detrended.csv"
YIELD_OBS_COL = 'yield_obs_detr'  # column in yield reference data

## Local folders with outputs
SIM_FOLDER = fr"D:\projects\KlimErtrag\out_remote_local\1_{RUN_ID}"
GRID_FOLDER = fr"D:\projects\KlimErtrag\raster\testing_with_ehsan\new_setup\{RUN_ID}"

## Output folder for figures
FIGURE_FOLDER = fr"D:\projects\KlimErtrag\out_remote_local\_figures\{RUN_ID}"

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

def find_event_in_lines_and_return_df(lines, event):

    start_ind = lines.index(f"{event}\n")
    col_ind = start_ind + 1
    end_ind = lines[start_ind:].index("\n") + start_ind
    lines = [item.strip('\n') for item in lines]
    lines = [item.split(',') for item in lines]
    values = lines[start_ind + 3:end_ind]
    df_dict = {col:[] for col in lines[col_ind]}
    for item in values:
        for i, col in enumerate(lines[col_ind]):
            df_dict[col].append(item[i])

    df = pd.DataFrame(df_dict)
    return df

def read_monica_results(pth, event):

    file = open(pth)
    lines = file.readlines()

    df = find_event_in_lines_and_return_df(lines=lines, event=event)
    df.drop(columns='CM-count', inplace=True)

    cols = df.columns.tolist()
    year_ind = cols.index('Year')
    del cols[year_ind]
    new_cols = ['Year'] + cols
    df = df[new_cols]

    return df

def merge_all_simulations(folder, event):
    pth_lst = glob.glob(f"{folder}/**/*.csv")

    df_lst = []
    for i, pth in enumerate(pth_lst):
        df = read_monica_results(pth=pth, event=event)
        if "Crop" in df.columns.tolist():
            df.drop(columns="Crop", inplace=True)

        df = df.T
        cols = df.iloc[0]
        values = df.iloc[1].tolist()
        df_lst.append(values)
        # df = df[1:]
        # df.columns = cols
        # df.reset_index(inplace=True, drop=True)
        # df.index = [i]
        # df = df.loc[~df.index.duplicated(keep='first')]
        # df_lst.append(df)

    cols = df.iloc[0]

    df_out = pd.DataFrame(df_lst)
    df_out.columns = cols

    # df_out = pd.concat(df_lst)
    df_out.index = range(len(df_out))
    df_out.dropna(inplace=True)
    for col in df_out.columns:
        df_out[col] = df_out[col].astype(float)
        df_out[col] = df_out[col].astype(int)

    print("Merging simulation dfs done!", folder)
    return df_out


def monica_grids_results_to_df(folder, event, mask_pth=None):

    event_to_var_name = {
        "Sowing": "sdoy",
        "Stage-2": "s2doy",
        "cereal-stem-elongation": "sedoy",
        "Stage-4": "s4doy",
        "Stage-5": "s5doy",
        "Stage-6": "s6doy",
        "Harvest": "hdoy"
    }

    var = event_to_var_name[event]
    pth_lst = glob.glob(f"{folder}\*{var}*.asc")

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
    print('Plotting done!', out_pth)

def merge_yearly_aggregated_simulated_yields_with_observations(df_sim):
    df_obs_ger = pd.read_csv(OBS_PTH_GER)
    df_obs_ger[YIELD_OBS_COL] = df_obs_ger[YIELD_OBS_COL] * 100

    year_lst = list(range(1999, 2020))
    min_year = min(year_lst)
    max_year = max(year_lst)
    df_obs_ger = df_obs_ger.loc[(df_obs_ger.Year >= min_year) & (df_obs_ger.Year <= max_year)]
    if CROP == "Winterroggen":
        crop = "Roggen"
    elif CROP == "Mais":
        crop = "Silomais"
    elif CROP == "Zucker-Ruebe":
        crop = "Zuckerrüben"
    else:
        crop = CROP
    df_obs_ger = df_obs_ger.loc[(df_obs_ger.Crop_name_de == crop)]

    df_obs_ger['Year'] = df_obs_ger['Year'].astype(int)
    df_ts_ger = pd.merge(df_obs_ger, df_sim, how='left', on=['Year'])

    return df_ts_ger

def plot_observat_simulation_time_series(df, obs_col, sim_col, min_year, max_year, title, out_pth):
    font_size = 12
    plt.rcParams['legend.title_fontsize'] = f'{font_size}'
    plt.rcParams["font.family"] = "Calibri"

    df = df.loc[(df['Year'] >= min_year) & df['Year'] <= max_year].copy()

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=cm2inch(16, 8))

    ax.plot(df['Year'], df[obs_col], c='red')
    ax.plot(df['Year'], df[sim_col], c='blue')

    ax.set_title(title, loc='left', fontsize=font_size)
    # ax.set_xlabel('Year')
    # ax.set_ylabel('Yield dry matter [kg/ha]')
    ax.set_xlabel('Jahr')
    ax.set_ylabel('Ertrag Trockenmasse [kg/ha]')

    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.01,
                     box.width, box.height * 0.95])

    custom_lines = [Line2D([0], [0], color='red', lw=1), Line2D([0], [0], color='blue', lw=1)]
    # legend_labels = ['Reference', 'Simulation']
    legend_labels = ['Referenz', 'Simulation']
    ax.legend(custom_lines, legend_labels, loc='upper center', bbox_to_anchor=(0.15, -0.08),
              ncol=2)

    plt.tight_layout()
    # plt.show()
    plt.savefig(out_pth)
    plt.close()

def plot_obs_vs_sim(df, obs_col, sim_col, title, out_pth):
    font_size = 12
    plt.rcParams['legend.title_fontsize'] = f'{font_size}'
    plt.rcParams["font.family"] = "Calibri"

    df = df[df[obs_col].notna()]
    df = df[df[sim_col].notna()]

    ## Calculate performance measures
    agr_ind = round(d_modified(df[obs_col], df[sim_col]), 2)
    rmse = round(math.sqrt(metrics.mean_squared_error(df[obs_col], df[sim_col])), 2)
    nmbe = round(np.mean(np.array(df[obs_col]) - np.array(df[sim_col])) / np.mean(np.array(df[obs_col])), 2)
    corr_matrix = np.corrcoef(np.array(df[obs_col]), np.array(df[sim_col]))
    corr = corr_matrix[0, 1]
    rsq = round(corr ** 2, 2)

    ## Plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=cm2inch(8, 8))

    ax.scatter(x=df[obs_col], y=df[sim_col], s=3)
    ax.annotate(text=f"nMBE: {nmbe}", xy=(0.051, 0.60), xycoords='axes fraction')
    ax.annotate(text=f"R²: {rsq}", xy=(0.051, 0.70), xycoords='axes fraction')
    ax.annotate(text=f"RMSE: {rmse}", xy=(0.051, 0.8), xycoords='axes fraction')
    ax.annotate(text=f"d_mod: {agr_ind}", xy=(0.051, 0.90), xycoords='axes fraction')

    max_val = max(df[sim_col].max(), df[obs_col].max())
    min_val = min(df[sim_col].min(), df[obs_col].min())

    ax.plot([min_val, max_val], [min_val, max_val], linewidth=0.5, c='black', linestyle='--')

    ax.set_title(title, loc='left', fontsize=font_size)
    # ax.set_xlabel('Reference')
    # ax.set_ylabel('Prediction')
    ax.set_xlabel('Referenz')
    ax.set_ylabel('Simulation')

    plt.tight_layout()
    # plt.show()
    plt.savefig(out_pth)
    plt.close()

def main():
    ###############################################
    ############ Phenology calibration ############
    ###############################################

    for event in EVENTS:
        try:
            df_sim = merge_all_simulations(folder=SIM_FOLDER, event=event)
            df_obs = pd.read_csv(OBS_PTH)

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

            phase_id = monica_to_dwd_stages[CROP][event]
            df_obs = df_obs.loc[df_obs['Phase_id'] == phase_id].copy()

            df_sim = pd.melt(df_sim, value_vars=df_sim.columns.tolist(), var_name='Year', value_name='value')
            df_sim['type'] = 'sim_sample'
            df_sim = df_sim[['Year', 'type', 'value']]

            df_obs = pd.melt(df_obs[['Referenzjahr', 'Jultag']], id_vars='Referenzjahr', value_vars='Jultag',
                             var_name='type', value_name='value')
            df_obs.columns = ['Year', 'type', 'value']
            df_obs['type'] = 'obs'

            df_plt = pd.concat([df_sim, df_obs])
            df_plt['Year'] = df_plt['Year'].astype(int)

            create_folder(FIGURE_FOLDER)
            out_pth = rf"{FIGURE_FOLDER}\{CROP}_{event}_boxplot.png"
            plot_simulations_vs_observations(df_plt, out_pth)
        except:
            pass

        # ### Test against grid simulations (comparison to Ehsans plots)
        # folder = GRID_FOLDER
        # if os.path.exists(folder):
        #     df_sim2 = monica_grids_results_to_df(folder, event)
        #     df_sim2 = pd.melt(df_sim2, value_vars=df_sim2.columns.tolist(), var_name='Year', value_name='value')
        #     df_sim2['type'] = 'sim_grid'
        #     df_sim2 = df_sim2[['Year', 'type', 'value']]
        #
        #     mask_pth = r"D:\projects\KlimErtrag\klimertrag\data\germany\dwd-stations-pheno_1000_25832_etrs89-utm32n.asc"
        #     df_sim3 = monica_grids_results_to_df(folder, event, mask_pth=mask_pth)
        #     df_sim3 = pd.melt(df_sim3, value_vars=df_sim3.columns.tolist(), var_name='Year', value_name='value')
        #     df_sim3['type'] = 'sim_grid_sample'
        #     df_sim3 = df_sim3[['Year', 'type', 'value']]
        #
        #     df_plt = pd.concat([df_sim, df_sim3, df_sim2, df_obs])
        #     df_plt['Year'] = df_plt['Year'].astype(int)
        #
        #     out_pth = rf"{FIGURE_FOLDER}\{CROP}_{event}_grid_simulations_boxplot_WB_dwd_pheno.png"
        #     plot_simulations_vs_observations(df_plt, out_pth)
        #
        #     ### Test against grid simulations (comparison to Ehsans plots)
        #     folder = GRID_FOLDER
        #     crop_mask_pth = r"D:\projects\KlimErtrag\raster\crop_maps\crop_specific_masks\CM_2017-2019_WB_1000m_25832_q3.asc"
        #     if os.path.exists(folder):
        #         df_sim3 = monica_grids_results_to_df(folder, event, crop_mask_pth)
        #         df_sim3 = pd.melt(df_sim3, value_vars=df_sim3.columns.tolist(), var_name='Year', value_name='value')
        #         df_sim3['type'] = 'sim_grid_masked'
        #         df_sim3 = df_sim3[['Year', 'type', 'value']]
        #         df_plt = pd.concat([df_sim, df_sim2, df_sim3, df_obs])
        #         df_plt['Year'] = df_plt['Year'].astype(int)
        #
        #         out_pth = rf"{FIGURE_FOLDER}\{CROP}_{event}_grid_simulations_boxplot.png"
        #         plot_simulations_vs_observations(df_plt, out_pth)

    ###############################################
    ############ FCM calibration ##################
    ###############################################
    create_folder(FIGURE_FOLDER)
    event = 'crop'
    df_sim = merge_all_simulations(folder=SIM_FOLDER, event=event)
    df_sim = df_sim.mean(axis=0)
    df_sim = pd.DataFrame(df_sim)
    df_sim['Year'] = df_sim.index
    df_sim.index = range(len(df_sim))
    df_sim.drop(index=21, inplace=True)
    df_sim.columns = ['simulated_value', 'Year']
    df_sim['Year'] = df_sim['Year'].astype(int)

    df_plt = merge_yearly_aggregated_simulated_yields_with_observations(df_sim)

    out_pth = rf"{FIGURE_FOLDER}\{CROP}_{event}_time_series_yields_1999-2019.png"
    plot_observat_simulation_time_series(df=df_plt,
                                             obs_col=YIELD_OBS_COL,
                                             sim_col='simulated_value',
                                             min_year=1999,
                                             max_year=2019,
                                             title=f'{CROP} Ø {1999}-{2019}',
                                             out_pth=out_pth)


    out_pth = rf"{FIGURE_FOLDER}\{CROP}_{event}_scatter-temporal_obs-sim_1999-2019.png"
    plot_obs_vs_sim(df=df_plt,
                    obs_col=YIELD_OBS_COL,
                    sim_col='simulated_value',
                    title=f'{CROP} Ø {1999}-{2019}',
                    out_pth=out_pth)


if __name__ == '__main__':
    main()