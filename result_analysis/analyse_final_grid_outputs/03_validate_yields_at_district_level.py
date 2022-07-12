import os
import numpy as np
from numpy.lib.function_base import diff
from numpy.lib.twodim_base import mask_indices
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D
import geopandas as gpd
from scipy.sparse import dok
from sklearn import metrics
import math
import glob

WD = "/beegfs/jaenicke/klimertrag/"
os.chdir(WD)

LDK_MSK = "/beegfs/jaenicke/klimertrag/result_analysis/data/administrative_masks/landkreise_1000.asc"  
LDK_SHP = "/beegfs/jaenicke/klimertrag/result_analysis/data/vector/GER_landkreise_25832.shp"  
CROP_MSK_PTH2 = "/beegfs/jaenicke/klimertrag/result_analysis/data/crop_masks/CM_2017-2019_{0}_1000m_25832_q3.asc"
# CROP_MSK_PTH2 = "/beegfs/jaenicke/klimertrag/result_analysis/data/crop_masks/CTM_17-19_mask_1000m_25832.asc"

YIELD_PTH = "/beegfs/jaenicke/klimertrag_temp/raster/sim-yields/11_final_results"

RUN_IDS = list(range(171,175))  #18,36,54,72,90,108,126,

OBS_PTH_LDK = f"/beegfs/jaenicke/klimertrag/result_analysis/data/reference_tables/yields/Ertraege_Landkreise_10_Frucharten_1999-2020_detrended.csv"
OBS_PTH_GER = f"/beegfs/jaenicke/klimertrag/result_analysis/data/reference_tables/yields/Ertraege_Deutschland_10_Frucharten_1999-2020_detrended.csv"

OUT_FOLDER_TABLES = "/beegfs/jaenicke/klimertrag_temp/tables/11_final_results"
OUT_FOLDER_FIGURES = "/beegfs/jaenicke/klimertrag_temp/figures/11_final_results/yield_validation"

YIELD_OBS_COL = 'yield_obs_detr'

## [Range of years to validation, Conversion factor fresh matter to dry matter, absolute max error for plotting]
CROPS = {
    'WW': [list(range(1999, 2020)),86],
    'SM': [list(range(1999, 2020)),32],
    'WB': [list(range(1999, 2020)),86, 2500],
    'WR': [list(range(1999, 2020)),91],
    'WRye': [list(range(1999, 2020)),86],
    'PO': [list(range(1999, 2020)),22.5],
    'SB': [list(range(1999, 2020)),86],
    'SU': [list(range(1999, 2020)),22.5]
}

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

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

def calc_perf_meas(df):
    print('')

    ##lst = [["MAE", "nMAE", "RMSE", "MBE", "nMBE", "d_modified", "R_squared"]]

    pred = 'yield_sim'
    ref = YIELD_OBS_COL #'yield_obs_detr'
    df = df.dropna(subset=[pred, ref])
    pred_data = df[pred].to_numpy()
    ref_data = df[ref].to_numpy()

    if len(pred_data) > 0:
        mae = round(metrics.mean_absolute_error(ref_data, pred_data), 3)
        nmae = round(metrics.mean_absolute_error(ref_data, pred_data) / np.mean(ref_data), 3)
        rmse = round(math.sqrt(metrics.mean_squared_error(ref_data, pred_data)), 3)
        mbe = np.mean(pred_data - ref_data)
        nmbe = np.mean(pred_data - ref_data) / np.mean(ref_data)
        agr_ind = round(d_modified(ref_data, pred_data), 3)

        corr_matrix = np.corrcoef(np.array(ref_data), np.array(pred_data))
        corr = corr_matrix[0, 1]
        rsq = round(corr ** 2, 2)

        out_lst = [mae, nmae, rmse, mbe, nmbe, agr_ind, rsq]

        return out_lst

    else:
        print('No predicted values in df.')

def plot_yield_maps(shp, obs_col, sim_col, title, out_pth):
    font_size = 12
    plt.rcParams['legend.title_fontsize'] = f'{font_size}'
    plt.rcParams["font.family"] = "DejaVu Sans"

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=cm2inch(20.0, 10))

    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right", size="5%", pad=0.1)

    max_val = max(shp[obs_col].max(), shp[sim_col].max())

    shp.plot(column=obs_col, ax=axs[0], vmin=0, vmax=max_val, cax=cax)
    shp.plot(edgecolor='black', facecolor="none", ax=axs[0], lw=0.1)
    shp.plot(column=sim_col, ax=axs[1], vmin=0, vmax=max_val, legend=True, cax=cax,
            #  legend_kwds={'label': "Yield dry matter [kg/ha]"})
            legend_kwds={'label': "Ertrag Trockenm. [kg/ha]"})
    shp.plot(edgecolor='black', facecolor="none", ax=axs[1], lw=0.1)

    # axs[0].set_title('Observed', loc='left')
    # axs[1].set_title('Simulated', loc='left')
    axs[0].set_title('Referenz', loc='left')
    axs[1].set_title('Simulation', loc='left')

    axs[0].set_axis_off()
    axs[1].set_axis_off()

    fig.suptitle(title)

    # plt.show()
    plt.tight_layout()
    plt.savefig(out_pth)
    plt.close()

def plot_yield_difference(shp, diff_col, title, out_pth, max_val=None):
    font_size = 12
    plt.rcParams['legend.title_fontsize'] = f'{font_size}'
    plt.rcParams["font.family"] = "DejaVu Sans"

    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=cm2inch(10.0, 11))

    divider = make_axes_locatable(axs)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    if not max_val:
        max_val = max(shp[diff_col].max(), abs(shp[diff_col].min()))
    else:
        max_val = max_val

    # shp.plot(column=diff_col, ax=axs, legend=True)
    shp.plot(column=diff_col, ax=axs, vmin=-max_val, vmax=max_val, legend=True, cax=cax,
            #  legend_kwds={'label': "Yield diff. dry matter [kg/ha]"}, cmap='RdYlBu')
            legend_kwds={'label': "Ertrag Trockenm [kg/ha]"}, cmap='RdYlBu')
    
    shp.plot(edgecolor='black', facecolor="none", ax=axs, lw=0.1)

    axs.set_title(title)

    axs.set_axis_off()

    # fig.suptitle(title)

    # plt.show()
    plt.tight_layout()
    plt.savefig(out_pth)
    plt.close()

def plot_obs_vs_sim(df, obs_col, sim_col, title, out_pth):
    font_size = 12
    plt.rcParams['legend.title_fontsize'] = f'{font_size}'
    plt.rcParams["font.family"] = "DejaVu Sans"

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

def plot_observat_simulation_time_series(df, obs_col, sim_col, min_year, max_year, title, out_pth):
    font_size = 12
    plt.rcParams['legend.title_fontsize'] = f'{font_size}'
    plt.rcParams["font.family"] = "DejaVu Sans"

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

def merge_observation_simulation_years_aggregated():
    ## Loop run_ids
    for run_id in RUN_IDS:

        yield_pth = f"{YIELD_PTH}/{run_id}"

        if not os.path.exists(yield_pth):
            print(f"Path to yield simulations with run-ID {run_id} does not exist. \n {yield_pth}. \n")
            return

        ## Get crop name of current run ID
        crop_name_dict = {
            'wheatwinterwheat':'WW',
            'barleywinterbarley':'WB',
            'ryewinterrye':'WRye',
            'maizesilagemaize':'SM',
            'sugarbeet':'SU',
            'potatomoderatelyearlypotato':'PO',
            'rapewinterrape':'WR',
            'barleyspringbarley':'SB'
        }
        search_term = rf'{yield_pth}/*Yield*.asc'
        file_lst = glob.glob(search_term)
        if file_lst:
            file_pth = file_lst[0]
            crop_name = os.path.basename(file_pth).split('_')[0]
            crop = crop_name_dict[crop_name]
        else:
            continue
        
        ## Get year list from global variables
        ## And frischmasse-trockenmasse conversion factor
        year_lst = CROPS[crop][0]
        conv_factor = CROPS[crop][1]
        yield_arr_lst = []

        ## Open Landkreise Raster and
        ldk_ras = gdal.Open(LDK_MSK)
        ldk_arr = ldk_ras.ReadAsArray()
        ndv_ldk = ldk_ras.GetRasterBand(1).GetNoDataValue()

        ldk_ids = np.unique(ldk_arr)
        ldk_ids = ldk_ids[ldk_ids != ndv_ldk]
        ldk_ids = ldk_ids[ldk_ids != 0]

        ## List for yields that are averaged over years
        sim_yields = []

        ## list for german wide means
        mean_ts1 = []
        ## list for mean of landkreise averages
        mean_ts2 = []
    
        for year in year_lst:
            print(year, crop)
            ## list all files of the current run
            search_term = rf'{yield_pth}/{crop_name}_Yield_{year}*.asc'
            file_pth = glob.glob(search_term)
            if file_pth:
                file_pth = file_pth[0]
            else:
                continue

            print(file_pth)

            ## Open simulated crop yields
            yield_ras = gdal.Open(file_pth)
            yield_arr = yield_ras.ReadAsArray()
            ndv_yields = yield_ras.GetRasterBand(1).GetNoDataValue()

            ## Create a mask that sets all no data values to 0
            ndv_mask = np.where(yield_arr == ndv_yields, 0, 1)

            ## Open crop mask aggregated
            crop_ras = gdal.Open(CROP_MSK_PTH2.format(crop))
            crop_arr = crop_ras.ReadAsArray()

            ## Create crop mask from crop array and ndv array
            mask = np.where(ndv_mask == 1, crop_arr, 0)

            ## Mask yield array with crop masked
            yield_arr_m = np.ma.masked_where(mask == 0, yield_arr)

            ## Convert dry matter to fresh matter
            ## !! Not necessary anymore, as reference data were already converted from fresh matter to dry matter !!
            # yield_arr_m = (yield_arr_m/conv_factor)*100

            yield_arr_lst.append(yield_arr_m.copy())

            ## Save mean german yield of current crop and year to list
            mean_yield = np.ma.mean(yield_arr_m)
            mean_ts1.append([crop, year, mean_yield])

        if not yield_arr_lst:
            print(f"There are no rasters for {crop}. Run ID: {run_id}")
            continue

        ## Average yields over all years
        yields_aggr = np.ma.mean(yield_arr_lst, axis=0)

        for fid in ldk_ids:
            ## Create landkreis-specific mask
            mask = np.where(ldk_arr == fid, 1, 0)

            ## Extract average yield for each landkreis across all years
            yields_aggr_m = np.ma.masked_where(mask == 0, yields_aggr)
            yields_aggr_m = np.ma.compressed(yields_aggr_m)
            sim_yield = np.mean(yields_aggr_m)
            sim_yields.append([fid, crop, sim_yield])

            ## Extract average yield for each landkreis and year separately
            for y, yield_arr in enumerate(yield_arr_lst):
                yield_arr_m = np.ma.masked_where(mask == 0, yield_arr)
                yield_arr_m = np.ma.compressed(yield_arr_m)
                mean_yield = np.ma.mean(yield_arr_m)
                year = year_lst[y]
                mean_ts2.append([crop, year, fid, mean_yield])

        ## Create df with mean yields of Germany
        cols = ['Crop', 'Year', 'mean_yield']
        df_sim_ger = pd.DataFrame(mean_ts1)
        df_sim_ger.columns = cols

        ## Create df with mean yields of landkreise averaged over all years
        cols = ['ID', 'Crop', 'yield_sim']
        df_sim = pd.DataFrame(sim_yields)
        df_sim.columns = cols

        ## Create df with mean yields of landkreises for each year separetely
        cols = ['Crop', 'Year', 'ID', 'mean_yield']
        df_sim_ldk = pd.DataFrame(mean_ts2)
        df_sim_ldk.columns = cols

        ## Open observed yields landkreis and for germany
        df_obs_l = pd.read_csv(OBS_PTH_LDK)
        df_obs_l[YIELD_OBS_COL] = df_obs_l[YIELD_OBS_COL] * 100
        df_obs_ger = pd.read_csv(OBS_PTH_GER)
        df_obs_ger[YIELD_OBS_COL] = df_obs_ger[YIELD_OBS_COL] * 100

        ## Change dtype of merge columns, so that merging works
        df_obs_l[YIELD_OBS_COL] = df_obs_l[YIELD_OBS_COL].apply(float)
        df_obs_l.Crop = df_obs_l.Crop.apply(str)
        df_obs_l.ID = df_obs_l.ID.apply(int)
        df_sim.Crop = df_sim.Crop.apply(str)
        df_sim.ID = df_sim.ID.apply(int)

        ## Subset reference dfs to period of interest
        min_year = min(year_lst)
        max_year = max(year_lst)
        df_obs_l = df_obs_l.loc[(df_obs_l.Year >= min_year) & (df_obs_l.Year <= max_year)]
        df_obs_ger = df_obs_ger.loc[(df_obs_ger.Year >= min_year) & (df_obs_ger.Year <= max_year)]

        #### For validation on Landkreise-level
        ## Average the yield over all years per landkreis
        df_obs_aggr = df_obs_l.groupby(['ID', 'Crop'])[YIELD_OBS_COL].mean().reset_index()

        ## Merge simulated with observed values
        df_val = pd.merge(df_sim, df_obs_aggr, how='left', on=['ID', 'Crop'])
        out_pth = rf"{OUT_FOLDER_TABLES}/{crop}_nuts3_GER_sim-obs_yields_years_aggr_{run_id}.csv"
        df_val.to_csv(out_pth, index=False)

        #### For validation on German-wide-level
        ## Average the yields over all landkreise per year just for comparison 
        ## (df with time series of mean mean yields of landkreise averages)
        df_sim_ger2 = df_sim_ldk.groupby(['Crop', 'Year'])['mean_yield'].mean().reset_index()

        ## Combine simulated yearly mean yields per (german-wide) with yearly mean yields of landkreise averages and
        ## observed yearly mean yields (landkreise averages)
        df_obs_ger['Year'] = df_obs_ger['Year'].astype(int)
        df_ts_ger = pd.merge(df_obs_ger, df_sim_ger, how='left', on=['Crop', 'Year'])
        df_ts_ger = pd.merge(df_ts_ger, df_sim_ger2, how='left', on=['Crop', 'Year'])
        df_ts_ger.columns = ['Crop', 'Crop_name_de', 'Year', 'yield_obs', 'yield_obs_detr', 'yield_sim_geravg', 'yield_sim_ldkavg']
        out_pth = rf"{OUT_FOLDER_TABLES}/{crop}_GER_sim-obs_yields_years_ts_{run_id}.csv"
        df_ts_ger.to_csv(out_pth, index=False)

        ## Combine time series of observed and simulated mean yields per landkreis
        df_obs_l['Year'] = df_obs_l['Year'].astype(int)
        df_ts_ldk = pd.merge(df_obs_l, df_sim_ldk, how='left', on=['Crop', 'Year', 'ID'])
        out_pth = rf"{OUT_FOLDER_TABLES}/{crop}_nuts3_GER_sim_vs_obs_yields_years_ts_{run_id}.csv"
        df_ts_ldk.to_csv(out_pth, index=False)

# def validate_results_year_specific():
#     df = pd.read_csv(r"tables\nuts3_GER_simulated_vs_observed_yields.csv")
#     shp = gpd.read_file(r"administrative\GER_landkreise_25832.shp")
#     shp.cca_2 = shp.cca_2.apply(int)

#     out_lst = []
#     for crop in CROPS:
#         year_lst = CROPS[crop][0]
#         for year in year_lst:
#             df_sub = df.loc[(df['crop'] == crop) & (df['year'] == year)]

#             ## Merge yields with shp
#             shp_plt = pd.merge(shp, df_sub, how='left', left_on='cca_2', right_on='fid')

#             ## Plot maps
#             out_pth = fr"figures\klimertrag_landkreise\comp_map_obs-sim_{crop}_{year}.png"
#             title=''
#             plot_yield_maps(shp_plt, obs_col=YIELD_OBS_COL, sim_col='yield_sim', title=title, out_pth=out_pth)

#             ## Validate accuracy
#             perf_meas = calc_perf_meas(df_sub)
#             out_lst.append([crop, year] + perf_meas)

#             ## Plot Scatterplot
#             out_pth = fr"figures\klimertrag_landkreise\scatter_obs-sim_{crop}_{year}.png"
#             plot_obs_vs_sim(df=df_sub, obs_col=YIELD_OBS_COL, sim_col='yield_sim', title=f'{crop} {year}', out_pth=out_pth)

#     cols = ["Crop", "Year", "MAE", "nMAE", "RMSE", "MBE", "nMBE", "d_modified", "R_squared"]
#     df_out = pd.DataFrame(out_lst)
#     df_out.columns = cols
#     out_pth = r"tables\validation_runs.csv"
#     df_out.to_csv(out_pth, index=False)


def validate_results_years_aggregated(run_id):

    yield_pth = f"{YIELD_PTH}/{run_id}"

    if not os.path.exists(yield_pth):
        print(f"Path to yield simulations with run-ID {run_id} does not exist. \n {yield_pth}. \n")
        return

    ## Get crop name of current run ID
    crop_name_dict = {
        'wheatwinterwheat':['WW','Winter-Weizen'],
        'barleywinterbarley':['WB','Winter-Gerste'],
        'ryewinterrye':['WRye','Winter-Roggen'],
        'maizesilagemaize':['SM','Silo-Mais'],
        'sugarbeet':['SU','Zuckerrübe'],
        'potatomoderatelyearlypotato':['PO','Kartoffel'],
        'rapewinterrape':['WR','Winter-Raps'],
        'barleyspringbarley':['SB','Sommer-Gerste']
    }
    search_term = rf'{yield_pth}/*.asc'
    file_lst = glob.glob(search_term)
    if file_lst:
        file_pth = file_lst[0]
        crop_name = os.path.basename(file_pth).split('_')[0]
        crop = crop_name_dict[crop_name][0]
        crop_ger = crop_name_dict[crop_name][1]
    else:
        return

    df = pd.read_csv(f"{OUT_FOLDER_TABLES}/{crop}_nuts3_GER_sim-obs_yields_years_aggr_{run_id}.csv")
    df_ts = pd.read_csv(f"{OUT_FOLDER_TABLES}/{crop}_GER_sim-obs_yields_years_ts_{run_id}.csv")
    shp = gpd.read_file(LDK_SHP)
    shp.cca_2 = shp.cca_2.apply(int)

    out_lst = []
    year_lst = CROPS[crop][0]
    min_year = min(year_lst)
    max_year = max(year_lst)

    ## ToDo: Move this into calculation of tables
    ## Calculate the difference between the observed and the simulated values
    df['diff_obsim'] = df[YIELD_OBS_COL] - df['yield_sim']

    df_sub = df.loc[(df['Crop'] == crop)]
    df_ts_sub = df_ts.loc[(df_ts['Crop'] == crop)]

    ## Merge yields with shp
    shp_plt = pd.merge(shp, df_sub, how='left', left_on='cca_2', right_on='ID')
    
    ## Remove Oberallgäu
    shp_plt.loc[shp_plt['cca_2'] == 9780, 'diff_obsim'] = None
    shp_plt.loc[shp_plt['cca_2'] == 9780, YIELD_OBS_COL] = None
    shp_plt.loc[shp_plt['cca_2'] == 9780, 'yield_sim'] = None

    ## Plot maps
    out_pth = fr"{OUT_FOLDER_FIGURES}/comp_map_obs-sim_{crop}_{min_year}-{max_year}_{run_id}.png"
    plot_yield_maps(shp_plt,
                    obs_col=YIELD_OBS_COL,
                    sim_col='yield_sim',
                    # title=f'{crop_ger} mean yields {min_year}-{max_year}',
                    title=f'{crop_ger} Ø {min_year}-{max_year}',
                    out_pth=out_pth)

    # max_val = CROPS[crop][2]
    out_pth = fr"{OUT_FOLDER_FIGURES}/diff_map_obs-sim_{crop}_{min_year}-{max_year}_{run_id}.png"
    plot_yield_difference(shp_plt,
                        diff_col='diff_obsim',
                        # title=f'{crop} - difference obs-sim of mean yields {min_year}-{max_year}',
                        title=f'{crop_ger} Ref-Sim Ø {min_year}-{max_year}',
                        out_pth=out_pth ) #  , max_val=max_val)

    ## Plot Scatterplot from spatial validation
    try:
        
        out_pth = fr"{OUT_FOLDER_FIGURES}/scatter-spatial_obs-sim_{crop}_{min_year}-{max_year}_{run_id}.png"
        plot_obs_vs_sim(df=df_sub,
                        obs_col=YIELD_OBS_COL,
                        sim_col='yield_sim',
                        # title=f'{crop} mean yields {min_year}-{max_year}',
                        title=f'{crop_ger} Ø {min_year}-{max_year}',
                        out_pth=out_pth)
    except:
        print(f"Probably no predicted or observed values for this time period. {min_year} - {max_year}.")

    ## Plot time series
    try:
        out_pth = rf"{OUT_FOLDER_FIGURES}/timeseries_obs-sim_{crop}_{min_year}-{max_year}_{run_id}.png"
        plot_observat_simulation_time_series(df=df_ts_sub,
                                                obs_col=YIELD_OBS_COL,
                                                sim_col='yield_sim_geravg',
                                                min_year=min_year,
                                                max_year=max_year,
                                                # title=f'{crop} mean yields {min_year}-{max_year}',
                                                title=f'{crop_ger} Ø {min_year}-{max_year}',
                                                out_pth=out_pth)
    except:
        print(f"Probably no predicted or observed values for this time period. {min_year} - {max_year}.")

    ## Plot Scatterplot from temporal validation
    try:
        out_pth = fr"{OUT_FOLDER_FIGURES}/scatter-temporal_obs-sim_{crop}_{min_year}-{max_year}_{run_id}.png"
        plot_obs_vs_sim(df=df_ts_sub,
                        obs_col=YIELD_OBS_COL,
                        sim_col='yield_sim_geravg',
                        # title=f'{crop} mean yields {min_year}-{max_year}',
                        title=f'{crop_ger} Ø {min_year}-{max_year}',
                        out_pth=out_pth)
    except:
        print(f"Probably no predicted or observed values for this time period. {min_year} - {max_year}.")

    try:
        ## Validate accuracy
        perf_meas = calc_perf_meas(df_sub)
        out_lst.append([crop] + perf_meas)

        cols = ["Crop", "MAE", "nMAE", "RMSE", "MBE", "nMBE", "d_modified", "R_squared"]
        df_out = pd.DataFrame(out_lst)
        df_out.columns = cols
        out_pth = fr"{OUT_FOLDER_TABLES}/validation_runs_years_aggregated_{run_id}.csv"
        df_out.to_csv(out_pth, index=False)
    except:
        print(f"Probably no predicted or observed values for this time period. {min_year} - {max_year}.")

if __name__ == '__main__':
    merge_observation_simulation_years_aggregated()
    for run_id in RUN_IDS:
        print(run_id)
        validate_results_years_aggregated(run_id=run_id)
    