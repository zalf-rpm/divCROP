import os
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import joblib
import geopandas as gpd

## Working directory
WD = "D:/projects/KlimErtrag/raster/11_final_results/results"

## Folder for maps and figures
FIGURES_FOLDER = "D:/projects/KlimErtrag/figures/klimertrag_landkreise/11_final_results"

## Germany shapefile for background
GERMANY_SHAPE = r"D:/projects/KlimErtrag/klimertrag/result_analysis/data/shapefiles/GER_bundeslaender.shp"

## Crop masks for masking results
# CROP_MSK_PTH = r"D:/projects/KlimErtrag/klimertrag/result_analysis/data/crop_masks/CM_2017-2019_{0}_1000m_25832_q3.asc"
CROP_MSK_PTH = r"D:/projects/KlimErtrag/klimertrag/result_analysis/data/crop_masks/CTM_17-19_mask_1000m_25832.asc"

OUTPUT_VAR_DICT = {
    "Yield" : {"title": "Ertrag",
               "axs_label" : "Ertrag Trockenmasse [kg/ha]"},
    "tradefavg" : {"title": "Transpirationsdefizit",
                   "axs_label" : "Transpirationsdefizit"},
    "frostredavg" : {"title": "Frostreduktion",
                     "axs_label": "Frostreduktion"}
                   }

SETUPS = {
    "wheatwinterwheat_klim":{
        "crop" : "wheatwinterwheat",
        "crop_name_ger" : "Winter-Weizen",
        "abbr": "WW",
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2099},
        "change_range" : {"min": -2500, "max" : 2500},
        "yield_max": 9000
    },
    "barleywinterbarley":{
        "crop" : "barleywinterbarley",
        "crop_name_ger" : "Wintergerste",
        "abbr": "WB",
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2099},
        "change_range" : {"min": -2500, "max" : 2500},
        "yield_max": 9000
    },
    "rapewinterrape":{
        "crop" : "rapewinterrape",
        "crop_name_ger" : "Winterraps",
        "abbr": "WR",
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2099},
        "change_range" : {"min": -2500, "max" : 2500},
        "yield_max": 5000
    },
    "maizesilagemaize": {
        "crop": "maizesilagemaize",
        "crop_name_ger": "Silo-Mais",
        "abbr": "SM",
        "hist_per": {"start": 1971, "end": 2000},
        "futur_per": {"start": 2031, "end": 2060},
        "futur_per2": {"start": 2070, "end": 2098},
        "change_range": {"min": -4500, "max": 4500},
        "yield_max": 22500
    },
    "wheatwinterwheat_sust":{
        "crop" : "wheatwinterwheat",
        "crop_name_ger" : "Winterweizen",
        "abbr": "WW",
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2099},
        "change_range" : {"min": -2500, "max" : 2500},
        "yield_max": 9000
    },
    "ryewinterrye":{
        "crop" : "ryewinterrye",
        "crop_name_ger" : "Winterroggen",
        "abbr": "WRye",
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2099},
            "change_range" : {"min": -2500, "max" : 2500},
        "yield_max": 9000
    },

    "barleyspringbarley":{
        "crop" : "barleyspringbarley",
        "crop_name_ger" : "Sommergerste",
        "abbr": "SB",
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2098},
        "change_range" : {"min": -2500, "max" : 2500},
        "yield_max": 5000
    },
    "potatomoderatelyearlypotato":{
        "crop" : "potatomoderatelyearlypotato",
        "crop_name_ger" : "Kartoffel",
        "abbr": "PO",
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2098},
        "change_range" : {"min": -5000, "max" : 5000},
        "yield_max": 12500
    },
    "sugarbeet":{
        "crop" : "sugarbeet",
        "crop_name_ger" : "Zuckerrübe",
        "abbr": "SU",
        "hist_per" : {"start": 1971, "end": 2000},
        "futur_per" : {"start": 2031, "end": 2060},
        "futur_per2" : {"start": 2070, "end": 2098},
        "change_range" : {"min": -9000, "max" : 9000},
        "yield_max": 22500
    },
    }

# SETUPS = {
#     "barleywinterbarley":{
#         "crop" : "barleywinterbarley",
#         "crop_name_ger" : "Winter-Gerste",
#         "abbr": "WB",
#         "hist_per" : {"start": 1971, "end": 2000},
#         "futur_per" : {"start": 2031, "end": 2060},
#         "futur_per2" : {"start": 2070, "end": 2099},
#         "change_range" : {"min": -2500, "max" : 2500},
#         "yield_max": 9000
#     }
#     }

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

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def get_extent(path):
    """
    Extracts the corners of a raster
	:param path: Path to raster including filename.
	:return: Minimum X, Minimum Y, Maximum X, Maximum Y
    """

    ds = gdal.Open(path)
    gt = ds.GetGeoTransform()
    width = ds.RasterXSize
    height = ds.RasterYSize
    minx = gt[0]
    miny = gt[3] + width * gt[4] + height * gt[5]
    maxx = gt[0] + width * gt[1] + height * gt[2]
    maxy = gt[3]

    return minx, maxx, miny, maxy

def plot_array(in_arr, extent, ignore_value, shp, title, out_pth, cmap=None, v_min=None, v_max=None):
    # in_arr = np.where(in_arr == ignore_value, np.nan, in_arr)
    in_arr = np.ma.masked_where(in_arr == ignore_value, in_arr)

    font_size = 20
    # plt.rcParams['legend.title_fontsize'] = f'{font_size}'
    plt.rcParams["font.family"] = "Calibri"
    plt.rcParams["font.size"] = f'{font_size}'

    if v_max and v_min == None:
        v_min = np.ma.min(in_arr)
    elif v_min and v_max == None:
        v_max = np.ma.max(in_arr)
    elif v_min == None and v_max == None:
        max_val = np.ma.max(in_arr)
        min_val = np.ma.min(in_arr)
        min_val_abs = abs(min_val)
        abs_max = max(max_val, min_val_abs)
        if min_val < 0:
            v_min = abs_max * -1
        else:
            v_min = min_val
        v_max = abs_max

    if not cmap:
        cmap = plt.cm.RdYlGn

    fig, ax = plt.subplots(figsize=cm2inch(25.0, 30))
    shp.plot(edgecolor='black', facecolor='lightgrey', ax=ax, zorder=0)
    img = ax.imshow(in_arr, cmap=cmap, extent=extent)
    ax.set_title(title, loc='left', fontsize=font_size)
    img.set_clim(vmin=v_min, vmax=v_max)
    fig.colorbar(img, ax=ax)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(out_pth, dpi=600)
    plt.close()

def plot_change_maps(sname, setup, wd, rcp, period1, period2, output_var, out_folder, ignore_value, v_min=None, v_max=None):
    create_folder(out_folder)

    crop = setup["crop"]

    start_fut = setup[period2]["start"]
    end_fut = setup[period2]["end"]
    start_hist = setup[period1]["start"]
    end_hist = setup[period1]["end"]

    in_pth = f"{wd}/change_maps/{sname}/{crop}_{rcp}_{output_var}_diff_{start_fut}-{end_fut}_{start_hist}-{end_hist}.tiff"
    in_arr = gdal.Open(in_pth).ReadAsArray()

    mask_arr = gdal.Open(CROP_MSK_PTH.format(setup["abbr"])).ReadAsArray()
    in_arr = np.where(mask_arr == 1, in_arr, ignore_value)

    shp = gpd.read_file(GERMANY_SHAPE)

    title = os.path.basename(in_pth).split('.')[0]
    out_pth = f"{out_folder}/{title}.png"

    extent = get_extent(in_pth)
    cmap = plt.cm.RdYlGn
    plot_array(in_arr=in_arr, extent=extent, ignore_value=ignore_value, shp=shp, title=title, out_pth=out_pth, cmap=cmap,
               v_min=v_min, v_max=v_max)

def plot_cov_maps(sname, setup, wd, rcp, period, output_var, statistic, out_folder, ignore_value, v_min=None, v_max=None):
    create_folder(out_folder)

    crop = setup["crop"]

    start = setup[period]["start"]
    end = setup[period]["end"]

    #"period_statistics/wheatwinterwheat_rcp26_Yield_1971-2000_cov.tiff"
    in_pth = f"{wd}/period_statistics/{sname}/{crop}_{rcp}_{output_var}_{start}-{end}_{statistic}.tiff"
    in_arr = gdal.Open(in_pth).ReadAsArray()

    mask_arr = gdal.Open(CROP_MSK_PTH.format(setup["abbr"])).ReadAsArray()
    in_arr = np.where(mask_arr == 1, in_arr, ignore_value)

    shp = gpd.read_file(GERMANY_SHAPE)

    title = os.path.basename(in_pth).split('.')[0]
    out_pth = f"{out_folder}/{title}.png"

    extent = get_extent(in_pth)
    cmap = plt.cm.YlOrRd
    plot_array(in_arr=in_arr, extent=extent, ignore_value=ignore_value, shp=shp, title=title, out_pth=out_pth, cmap=cmap,
               v_min=v_min, v_max=v_max)

def plot_array_grid(sname, setup, wd, rcp, period, output_var, statistic, title_descr, out_folder, ignore_value,
                    cmap=None, v_min=None, v_max=None):

    create_folder(out_folder)

    crop = setup["crop"]
    crop_name = setup["crop_name_ger"]
    start = setup[period]["start"]
    end = setup[period]["end"]

    title = f'{crop_name} {OUTPUT_VAR_DICT[output_var]["title"]} {statistic}{title_descr}'
    out_pth = f"{out_folder}/{crop_name}_{rcp}_{OUTPUT_VAR_DICT[output_var]['title']}_{start}-{end}_{statistic}.png"

    year_lst = list(range(start, end+1))

    mask_arr = gdal.Open(CROP_MSK_PTH.format(setup["abbr"])).ReadAsArray()
    arr_lst = []
    for year in year_lst:
        #barleywinterbarley_rcp26_frostredavg_1971_mean.tiff
        in_pth = f"{wd}/yearly_averages/{sname}/{crop}_{rcp}_{output_var}_{year}_{statistic}.tiff"
        arr = gdal.Open(in_pth).ReadAsArray()
        arr = np.where(mask_arr == 1, arr, ignore_value)
        arr = np.ma.masked_where(arr == ignore_value, arr)
        arr_lst.append(arr)

    extent = get_extent(in_pth)

    ## determine vmin and vmax
    min_lst = [np.ma.min(arr) for arr in arr_lst]
    max_lst = [np.ma.max(arr) for arr in arr_lst]

    if v_max and v_min == None:
        v_min = min(min_lst)
    elif v_min and v_max == None:
        v_max = max(max_lst)
    elif v_min == None and v_max == None:
        max_val = max(max_lst)
        min_val = min(min_lst)
        min_val_abs = abs(min(min_lst))
        abs_max = max(max_val, min_val_abs)
        if min_val < 0:
            v_min = abs_max * -1
        else:
            v_min = min_val
        v_max = abs_max

    shp = gpd.read_file(GERMANY_SHAPE)

    font_size = 11
    # plt.rcParams['legend.title_fontsize'] = f'{font_size}'
    plt.rcParams["font.family"] = "Calibri"
    plt.rcParams["font.size"] = f'{font_size}'

    fig, axs = plt.subplots(nrows=3, ncols=10, figsize=cm2inch(90, 30))

    if not cmap:
        cmap = plt.cm.RdYlGn

    for i, arr in enumerate(arr_lst):

        ix = np.unravel_index(i, axs.shape)

        shp.plot(edgecolor='black', facecolor='lightgrey', ax=axs[ix], zorder=0)
        img = axs[ix].imshow(arr, cmap=cmap, extent=extent)
        axs[ix].set_title(year_lst[i], loc='left', fontsize=font_size)
        img.set_clim(vmin=v_min, vmax=v_max)
        fig.colorbar(img, ax=axs[ix])
        axs[ix].axis("off")

    # for i in range(axs.shape[0]):
    #     axs[i, 0].set_ylabel(OUTPUT_VAR_DICT[output_var]["axs_label"])

    fig.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_pth, dpi=600)
    plt.close()

def plot_period_comp(sname, setup, wd, rcp, output_var, statistic, title_descr, out_folder, ignore_value,
                     v_min=None, v_max=None):
    create_folder(out_folder)

    crop = setup["crop"]

    arr_lst = []
    label_lst = []
    for period in ["hist_per", "futur_per", "futur_per2"]:
        start = setup[period]["start"]
        end = setup[period]["end"]
        in_pth = f"{wd}/period_statistics/{sname}/{crop}_{rcp}_{output_var}_{start}-{end}_{statistic}.tiff"
        arr = gdal.Open(in_pth).ReadAsArray()
        mask_arr = gdal.Open(CROP_MSK_PTH.format(setup["abbr"])).ReadAsArray()
        arr = np.where(mask_arr == 1, arr, ignore_value)
        arr_lst.append(arr)
        label = f"{start}-{end}"
        label_lst.append(label)

    shp = gpd.read_file(GERMANY_SHAPE)

    title = f"{crop}_{rcp}_{output_var}_{statistic}_{title_descr}"
    out_pth = f"{out_folder}/{title}.png"

    extent = get_extent(in_pth)
    cmap = plt.cm.YlOrBr

    def plot_multiple_arrays(in_arr_lst, label_lst, extent, ignore_value, shp, title, out_pth, cmap=cmap, v_min=None,
                             v_max=None):
        font_size = 20
        # plt.rcParams['legend.title_fontsize'] = f'{font_size}'
        plt.rcParams["font.family"] = "Calibri"
        plt.rcParams["font.size"] = f'{font_size}'

        num_labels = len(label_lst)
        num_arr = len(in_arr_lst)
        if num_labels < num_arr:
            print("Not enough labels provided! Adding empty strings to label list.")
            for i in num_arr - num_labels:
                label_lst.append('')

        ## determine vmin and vmax
        arr_lst = [np.ma.masked_where(in_arr == ignore_value, in_arr) for in_arr in in_arr_lst]
        min_lst = [np.ma.min(arr) for arr in arr_lst]
        max_lst = [np.ma.max(arr) for arr in arr_lst]

        if v_max and v_min == None:
            v_min = min(min_lst)
        elif v_min and v_max == None:
            v_max = max(max_lst)
        elif v_min == None and v_max == None:
            max_val = max(max_lst)
            min_val = min(min_lst)
            min_val_abs = abs(min(min_lst))
            abs_max = max(max_val, min_val_abs)
            if min_val < 0:
                v_min = abs_max * -1
            else:
                v_min = min_val
            v_max = abs_max

        ## create plot
        ncols = len(in_arr_lst)

        fig, axs = plt.subplots(nrows=1, ncols=ncols, figsize=cm2inch(25.0 * ncols, 30))

        if not cmap:
            cmap = plt.cm.YlOrBr

        for i, arr in enumerate(arr_lst):
            ix = np.unravel_index(i, axs.shape)

            shp.plot(edgecolor='black', facecolor='lightgrey', ax=axs[ix], zorder=0)
            img = axs[ix].imshow(arr, cmap=cmap, extent=extent)
            axs[ix].set_title(label_lst[i], loc='left', fontsize=font_size)
            img.set_clim(vmin=v_min, vmax=v_max)
            fig.colorbar(img, ax=axs[ix])
            axs[ix].axis("off")

        fig.suptitle(title)
        plt.tight_layout()
        plt.savefig(out_pth, dpi=300)
        plt.close()
        print(out_pth, "done!")

    plot_multiple_arrays(in_arr_lst=arr_lst, label_lst=label_lst, extent=extent, ignore_value=ignore_value, shp=shp,
                         title=title, out_pth=out_pth, cmap=None, v_min=v_min, v_max=v_max)

def plot_boxplot_grid(wd, setups, output_var, statistic, ignore_value, title, out_folder):
    create_folder(out_folder)
    out_pth = f"{out_folder}/{output_var}_{statistic}_grid_rcp_period_comp.png"

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=cm2inch(30,30))

    for i, sname in enumerate(setups):

        setup = setups[sname]

        ix = np.unravel_index(i, axs.shape)

        crop = setup["crop"]
        crop_name = setup["crop_name_ger"]

        mask_arr = gdal.Open(CROP_MSK_PTH.format(setup["abbr"])).ReadAsArray()
        data_lst = []
        color_labels = []
        for period in ["hist_per", "futur_per", "futur_per2"]:
            arr_lst = []
            for rcp in ["rcp26", "rcp45", "rcp85"]:
                start = setup[period]["start"]
                end = setup[period]["end"]
                in_pth = f"{wd}/period_statistics/{sname}/{crop}_{rcp}_{output_var}_{start}-{end}_{statistic}.tiff"
                arr = gdal.Open(in_pth).ReadAsArray()
                arr = np.where(mask_arr == 1, arr, ignore_value)
                arr = arr.flatten()
                arr = arr[arr != ignore_value]
                arr_lst.append(arr)
            data_lst.append(arr_lst)
            label = f"{start}-{end}"
            color_labels.append(label)

        data_a = data_lst[0]
        data_b = data_lst[1]
        data_c = data_lst[2]

        ticks = ['rcp26', 'rcp45', 'rcp85']

        def set_box_color(bp, color):
            plt.setp(bp['boxes'], color=color)
            plt.setp(bp['whiskers'], color=color)
            plt.setp(bp['caps'], color=color)
            plt.setp(bp['medians'], color=color)

        bpl = axs[ix].boxplot(data_a, positions=np.array(range(len(data_a))) * 3.0 - 0.65, sym='', widths=0.6)
        bpm = axs[ix].boxplot(data_b, positions=np.array(range(len(data_b))) * 3.0 + 0.0, sym='', widths=0.6)
        bpr = axs[ix].boxplot(data_c, positions=np.array(range(len(data_b))) * 3.0 + 0.65, sym='', widths=0.6)
        set_box_color(bpl, '#D7191C')  # colors are from http://colorbrewer2.org/
        set_box_color(bpm, '#2C7BB6')
        set_box_color(bpr, 'black')

        # draw temporary red and blue lines and use them to create a legend
        axs[ix].plot([], c='#D7191C', label=color_labels[0])
        axs[ix].plot([], c='#2C7BB6', label=color_labels[1])
        axs[ix].plot([], c='black', label=color_labels[2])
        axs[ix].legend()

        axs[ix].set_xticks(range(0, len(ticks) * 3, 3))
        axs[ix].set_xticklabels(ticks)
        axs[ix].set_xlim(-3, len(ticks) * 3)
        axs[ix].set_title(crop_name)

    for i in range(axs.shape[0]):
        axs[i, 0].set_ylabel(OUTPUT_VAR_DICT[output_var]["axs_label"])
    for i in range(axs.shape[1]):
        axs[axs.shape[0]-1, i].set_xlabel("RCP-Szenario")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.suptitle(title)
    # plt.show()
    plt.savefig(out_pth)
    plt.close()

def plot_boxplot_comp(sname, setup, wd, output_var, statistic, title_descr, out_folder, ignore_value):
    #plot_period_comp(setup, wd, rcp, output_var, statistic, title_descr, out_folder, ignore_value, v_min=None, v_max=None)

    create_folder(out_folder)

    crop = setup["crop"]
    crop_name = setup["crop_name_ger"]

    mask_arr = gdal.Open(CROP_MSK_PTH.format(setup["abbr"])).ReadAsArray()
    data_lst = []
    label_lst = []
    for period in ["hist_per", "futur_per", "futur_per2"]:
        arr_lst = []
        for rcp in ["rcp26", "rcp45", "rcp85"]:
            start = setup[period]["start"]
            end = setup[period]["end"]
            in_pth = f"{wd}/period_statistics/{sname}/{crop}_{rcp}_{output_var}_{start}-{end}_{statistic}.tiff"
            arr = gdal.Open(in_pth).ReadAsArray()
            arr = np.where(mask_arr == 1, arr, ignore_value)
            arr = arr.flatten()
            arr = arr[arr != ignore_value]
            arr_lst.append(arr)
        data_lst.append(arr_lst)
        label = f"{start}-{end}"
        label_lst.append(label)

    title = f"{crop_name} {output_var} {statistic}{title_descr}"
    out_pth = f"{out_folder}/{crop}_{output_var}_{statistic}_rcp_period_comp.png"

    def plot_boxplots(data_a, data_b, data_c, ticks, color_labels, title, out_pth):
        # ticks = ['rcp26', 'rcp45', 'rcp85']
        # color_labels = ['1970-2000', '2031-2060', '2070-2099']

        def set_box_color(bp, color):
            plt.setp(bp['boxes'], color=color)
            plt.setp(bp['whiskers'], color=color)
            plt.setp(bp['caps'], color=color)
            plt.setp(bp['medians'], color=color)

        plt.figure()

        bpl = plt.boxplot(data_a, positions=np.array(range(len(data_a))) * 3.0 - 0.65, sym='', widths=0.6)
        bpm = plt.boxplot(data_b, positions=np.array(range(len(data_b))) * 3.0 + 0.0, sym='', widths=0.6)
        bpr = plt.boxplot(data_c, positions=np.array(range(len(data_b))) * 3.0 + 0.65, sym='', widths=0.6)
        set_box_color(bpl, '#D7191C')  # colors are from http://colorbrewer2.org/
        set_box_color(bpm, '#2C7BB6')
        set_box_color(bpr, 'black')

        # draw temporary red and blue lines and use them to create a legend
        plt.plot([], c='#D7191C', label=color_labels[0])
        plt.plot([], c='#2C7BB6', label=color_labels[1])
        plt.plot([], c='black', label=color_labels[2])
        plt.legend()

        plt.xticks(range(0, len(ticks) * 3, 3), ticks)
        plt.xlim(-3, len(ticks) * 3)
        # plt.ylim(0, 8)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.suptitle(title)
        # plt.show()
        plt.savefig(out_pth)
        plt.close()

    plot_boxplots(data_a=data_lst[0], data_b=data_lst[1], data_c=data_lst[2], ticks=['rcp26', 'rcp45', 'rcp85'],
                  color_labels=label_lst, title=title, out_pth=out_pth)

def work_func(sname):
    output_var = "Yield"

    print(sname, output_var)
    setup = SETUPS[sname]

    ################################################################################################################
    ############ Boxplots of yield - Single crops comparison of periods and rcp scenarios ##########################
    ################################################################################################################

    out_folder = f"{FIGURES_FOLDER}/boxplots/{sname}"
    plot_boxplot_comp(sname, setup, WD, output_var, "mean", '', out_folder, -9999)

    for rcp in ["rcp26", "rcp45", "rcp85"]:
        print(rcp)
        ############################################################################################################
        ############ Change maps - difference between two periods ##################################################
        ############################################################################################################

        out_folder = f"{FIGURES_FOLDER}/change_maps/{sname}/{output_var}"
        v_min = setup["change_range"]["min"]
        v_max = setup["change_range"]["max"]
        plot_change_maps(sname, setup, WD, rcp, "hist_per", "futur_per", output_var, out_folder, ignore_value=0,
                         v_min=v_min, v_max=v_max)
        plot_change_maps(sname, setup, WD, rcp, "hist_per", "futur_per2", output_var, out_folder, ignore_value=0,
                         v_min=v_min, v_max=v_max)

        ############################################################################################################
        ############ Comparison of yields across periods with different statistics #################################
        ############################################################################################################

        out_folder = f"{FIGURES_FOLDER}/period_statistics/{sname}/{output_var}_mean"
        v_max = setup["yield_max"]
        plot_period_comp(sname, setup, WD, rcp, output_var, "mean", "period_comp", out_folder, ignore_value=-9999,
                         v_min=0, v_max=v_max)
        out_folder = f"{FIGURES_FOLDER}/period_statistics/{sname}/{output_var}_std"
        plot_period_comp(sname, setup, WD, rcp, output_var, "std", "period_comp", out_folder, ignore_value=-9999,
                         v_min=None, v_max=None)
        out_folder = f"{FIGURES_FOLDER}/period_statistics/{sname}/{output_var}_cov"
        plot_period_comp(sname, setup, WD, rcp, output_var, "cov", "period_comp", out_folder, ignore_value=-9999,
                         v_min=None, v_max=None)

        for period in ["hist_per", "futur_per", "futur_per2"]:
            print(period)
            ########################################################################################################
            ############ One map per crop, rcp, period, variable and statistic #####################################
            ########################################################################################################

            out_folder = f"{FIGURES_FOLDER}/period_statistics/{sname}"
            plot_cov_maps(sname, setup, WD, rcp, period, output_var, "cov", out_folder, ignore_value=-9999,
            v_min=0, v_max=0.5)

            ########################################################################################################
            ############ Time series of maps per crop, rcp, period, variable and statistic #########################
            ########################################################################################################

            # out_folder = f"{FIGURES_FOLDER}/map_time_series/{sname}"
            # plot_array_grid(sname=sname, setup=setup, wd=WD, rcp=rcp, period=period, output_var=output_var,
            #                 statistic="mean", title_descr="", out_folder=out_folder, ignore_value=-9999,
            #                 cmap=None, v_min=None, v_max=None)


    # for output_var in ["tradefavg", "frostredavg"]:
    #     print(output_var)
    #
    #     print(sname, output_var)
    #     setup = SETUPS[sname]
    #
    #     ############################################################################################################
    #     ############ Boxplots of variable - Single crops comparison of periods and rcp scenarios ###################
    #     ############################################################################################################
    #
    #     ## plot boxplots
    #     # out_folder = f"{FIGURES_FOLDER}/boxplots/{sname}"
    #     # plot_boxplot_comp(sname, setup, WD, output_var, "mean", '', out_folder, -9999)
    #
    #     for rcp in ["rcp26", "rcp45", "rcp85"]:
    #         print(rcp)
    #
    #         ########################################################################################################
    #         ############ Change maps - difference between two periods ##############################################
    #         ########################################################################################################
    #         # out_folder = f"{FIGURES_FOLDER}/change_maps/{sname}/{output_var}"
    #         # plot_change_maps(sname, setup, WD, rcp, "hist_per", "futur_per", output_var, out_folder, ignore_value=0,
    #         #                  v_min=-1, v_max=1)
    #         # plot_change_maps(sname, setup, WD, rcp, "hist_per", "futur_per2", output_var, out_folder,
    #         #                  ignore_value=0, v_min=-1, v_max=1)
    #
    #         ########################################################################################################
    #         ############ Comparison of variable across periods with different statistics ###########################
    #         ########################################################################################################
    #
    #         out_folder = f"{FIGURES_FOLDER}/period_statistics/{sname}/{output_var}_mean"
    #         plot_period_comp(sname, setup, WD, rcp, output_var, "mean", "period_comp", out_folder,
    #                          ignore_value=-9999, v_min=None, v_max=None)
    #         out_folder = f"{FIGURES_FOLDER}/period_statistics/{sname}/{output_var}_std"
    #         plot_period_comp(sname, setup, WD, rcp, output_var, "std", "period_comp", out_folder,
    #                          ignore_value=-9999, v_min=None, v_max=None)
    #         out_folder = f"{FIGURES_FOLDER}/period_statistics/{sname}/{output_var}_cov"
    #         plot_period_comp(sname, setup, WD, rcp, output_var, "cov", "period_comp", out_folder,
    #                          ignore_value=-9999, v_min=None, v_max=None)
    #
    #         for period in ["hist_per", "futur_per", "futur_per2"]:
    #             print(period)
    #
    #             ####################################################################################################
    #             ############ One map per crop, rcp, period, variable and statistic #################################
    #             ####################################################################################################
    #
    #             # out_folder = f"{FIGURES_FOLDER}/period_statistics/{sname}"
    #             # plot_cov_maps(sname, setup, WD, rcp, period, output_var, "cov", out_folder, ignore_value=-9999,
    #             #               v_min=0, v_max=0.5)

def main():
    ####################################################################################################################
    ############ Boxplot grid of yield - Comparison of crops, periods and rcp scenarios ################################
    ####################################################################################################################
    for output_var in ["Yield"]:

        out_folder = f"{FIGURES_FOLDER}/boxplots"
        plot_boxplot_grid(wd=WD, setups=SETUPS, output_var=output_var, statistic="mean", ignore_value=-9999,
                          title=OUTPUT_VAR_DICT[output_var]["title"], out_folder=out_folder)

if __name__ == '__main__':
    main()

    # joblib.Parallel(n_jobs=1)(joblib.delayed(work_func)(sname) for sname in SETUPS)