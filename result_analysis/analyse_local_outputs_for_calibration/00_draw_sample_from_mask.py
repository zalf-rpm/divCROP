import numpy as np
from osgeo import gdal
import os

WD = r"D:\projects\KlimErtrag\klimertrag\result_analysis\data"

N = 10000

MASK_PTH = fr"{WD}\crop_masks\CM_2017-2019_WB_1000m_25832_q3.asc"
OUT_PTH = fr"{WD}\random_samples\CM_2017-2019_WB_1000m_25832_q3_random_sample_{N}.tiff"

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

def main():

    ## Open Mask
    ras = gdal.Open(MASK_PTH)
    arr = ras.ReadAsArray()
    gt = ras.GetGeoTransform()
    pr = ras.GetProjection()

    ## Get tuple of x and y indeces where mask is 1.
    ind_mask = np.where(arr == 1)

    ## Get a random sample in the range of the length of the tuple
    rand_num = np.random.choice(len(ind_mask[0]), N, replace=False).tolist()

    ## Get random indeces of mask
    rand_inds = list(zip(ind_mask[0][rand_num], ind_mask[1][rand_num]))

    ## Create new mask
    new_mask = np.full(arr.shape, 0)
    for coord in rand_inds:
        new_mask[coord] = 1

    write_array_to_raster(in_array=new_mask,
                          out_path=OUT_PTH,
                          gt=gt,
                          pr=pr,
                          no_data_value=-9999,
                          type_code=None,
                          options=['COMPRESS=DEFLATE', 'PREDICTOR=1'],
                          driver='GTiff')

    os.system(f"gdal_translate -of AAIGrid {OUT_PTH} {OUT_PTH[:-4]}.asc")
    os.remove(f'{OUT_PTH}')

if __name__ == '__main__':
    main()