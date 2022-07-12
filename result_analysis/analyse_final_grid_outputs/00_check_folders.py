import os
import glob
from osgeo import gdal

wd = "/beegfs/jaenicke/klimertrag/raster/sim-yields/10_predictions_future"
os.chdir(wd)

lst = glob.glob("/beegfs/jaenicke/klimertrag/raster/sim-yields/10_predictions_future/*")

for i in lst:
    s_term =  f"{i}/*.asc"
    s_lst = glob.glob(s_term)
    print(i, len(s_lst))

task_lst = [str(i) for i in range(1, 154)]
done_lst = [os.path.basename(pth) for pth in lst]
miss_lst = [int(i) for i in task_lst if i not in done_lst]
print(miss_lst)

####################################
########### Check files ############

for i in lst:
    s_term =  f"{i}/*.asc"
    s_lst = glob.glob(s_term)
    for pth in s_lst:
        try:
            arr = gdal.Open(pth).ReadAsArray()
        except:
            print(i, pth)