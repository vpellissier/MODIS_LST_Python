import tempfile
import pymodis
from osgeo import gdal
import os
from calendar import monthrange
import numpy as np
from shutil import rmtree
import csv
import urllib.request
import io
import datetime
from joblib import Parallel, delayed
import itertools
import glob

'''
mosaics_parameters
Return a list containing every unique combination of Product, Year, Month, day_night.
Each element is a list with a unique combination.
The lists starts at the earliest month for which HDF files are available through the whole month
(March 2000 for MOD11A2 and August 2002 for MYD11A2) and stops one month before the current month.
'''
def mosaics_parameters():
    prod = ["MOD11A2.006", "MYD11A2.006"]
    current_month = datetime.datetime.now().month
    current_year = datetime.datetime.now().year
    year = list(range(2000, current_year + 1))
    month = list(range(1,13))
    day_night = ['day', 'night']
    combinations = list(itertools.product(prod, year, month, day_night))
    #removing every combinations not existing (before March 2000 for MOD11A2
    # and before July 2002 for MYD11A2)
    combinations = [row for row in combinations 
    if not (row[1] == 2018 and row[2] >= current_month)]
    combinations = [row for row in combinations 
    if not (row[0] == 'MOD11A2.006' and row[1] == 2000 and row[2] < 3)]
    combinations = [row for row in combinations 
    if not (row[0] == 'MYD11A2.006' and row[1] < 2002)]
    combinations = [row for row in combinations 
    if not (row[0] == 'MYD11A2.006' and row[1] == 2002 and row[2] < 7)]
    return(combinations)

'''
existing_tiles
Return an array containing all the tiles present, under
the forme [[tilesH], [tilesV]], with tilesH the horizontal identifiers
of the tile and tilesV the vertical identifiers
'''
def existing_tiles():
    url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MOD11A2/2005/001.csv'
    webpage = urllib.request.urlopen(url)
    datareader = csv.reader(io.TextIOWrapper(webpage))
    data = list(datareader)
    tilesH = []
    tilesV = []
    for row in data[1:]:
        tilesH.append(int(row[0][18:20]))
        tilesV.append(int(row[0][21:23]))
    return([tilesH, tilesV])


''' 
save_raster
Takes a numpy array along with georeferencing informations (projection and
transformation) and save it as a geotiff in the specified directory

Parameters:
value: ndarray containing the values to be written as raster
transfo: transformation information for the raster
projection: projection of the raster
save_path: full path of the raster file ('path/raster_dir/raster.tif')
in_na_value: no_data value of the input array
out_na_value: no_data value of the output GeoTiff
'''
def save_raster(value, transfo, projection, 
save_path, in_na_value, out_na_value = 0):
    # Replacing the no_data value from the input data (might be NAs) by output no_data value (default 0)
    if in_na_value == 'NA':
        value[np.isnan(value)] = out_na_value
    else:
        value[value == in_na_value] = out_na_value
    driver = gdal.GetDriverByName('GTiff')
    nrows = value.shape[0]
    ncols = value.shape[1]
    # Opening a raster
    out_raster = driver.Create(save_path, ncols, nrows, 1, 
        gdal.GDT_UInt16, ['COMPRESS=LZW', 'TILED=YES'])
    # Setting georeferencing attributes
    out_raster.SetGeoTransform(transfo)
    out_raster.SetProjection(projection)
    # Filling the raster with value
    out_raster.GetRasterBand(1).WriteArray(value)
    # Setting the no_data value
    out_raster.GetRasterBand(1).SetNoDataValue(out_na_value)
    # Closing the raster to write it on the HD
    out_raster = None


'''
raster_high_qc
Takes two raster (well, numpy arrays) as input (LST and QC) and return an array with
only the highest quality pixels (QA = 00).

Parameters:
LST: the LST raster (ndarray)
QC: the Quality Control raster matching the LST raster (ndarray)
'''
def raster_high_qc(LST, QC):
    '''
    Creates a mask with the QC raster. Here binary with  
    bits 1 and 0 equal to 0 correspond to the integers 0, 4 and 8
    '''
    QC_logical = np.logical_or(QC == 0, QC == 4, QC == 8) # True if high QC pixel, False otherwise
    # An array with LST value for high quality pixels
    HQ_LST = QC_logical * LST 
    # An array indicating wether the data for this pixel is missing
    HQ_LST = HQ_LST.astype(int)
    non_missing = HQ_LST != 0 
    return([HQ_LST, non_missing])



'''
df_month_hdf
Download HDFs of a tile for a whole month in a year (either day or night and either Terra or Aqua)
and extracts two numpy arrays (corresponding to the LST and QC layer) as well as 
the associated georeferencing inportmations (projection and transformation).
Then, the arrays are processed by raster_high_qc so that only the pixel with good QA bits are kept and
averaged to create a monthly composite of the tile. 
The resulting array is saved as a temporary GeoTiff using the following naming
convention 'xxxx.Ayyyymm.tile.day_night.tif' (xxxx: product name, yyyy: year,
    mm: month', tile: tile ID, day_night: time of the day).
This function could be modified to implement a form of temporal gap filling by (1) 
downloading HDF before and after the beginning of the month, (2) creating cloud masks,
(3) for each array within the month, averaging the LST value of the cloudy pixel 
from the array located before and after.

Parameters:
product (str): Product to dowload, either 'MOD11A2.006' or 'MYD11A2.006'
day_night (str): Time of the day, either 'day' or 'night'
tileH (int): horizontal ID of the tile (from 0 to 35)
tileV (int): vertical ID of the tile (from = to 17)
year (int): year for which the data must be processed
month (int): month for which the composite must be compiled
usr (str): user name as entered for https://e4ftl01.cr.usgs.gov
pwd (str): password as entered for https://e4ftl01.cr.usgs.gov
temp_path_geotiff (str): root path in which the monthly composite for the tile is saved
    (ex: 'C:/Geotiff'). The resulting composites will be saved in a directory
    'save_path/xxx.Ayyyymm.day_night (xxxx: product name, yyyy: year,
    mm: month').
'''
def dl_month_hdf(product, day_night, tileH, tileV, 
    year, month, usr, pwd, temp_path_geotiff):
    # Defining product path based on product name
    if product == "MOD11A2.006":
        path = "MOLT"
    else:
        path = "MOLA"
    
    # Defining layers to extract (day vs night)
    if day_night == "day":
        lst_sds_id, qc_sds_id = (0,1)
    else:
        lst_sds_id, qc_sds_id = (4,5)
    # Defining first and last day of the month
    tile = 'h' + str(tileH).zfill(2) + 'v' + str(tileV).zfill(2)
    last_day = monthrange(year, month)[1]
    year = str(year)
    month = str(month).zfill(2)
    first_date = year + '-' + month + '-' + '01'
    last_date = year + '-' + month + '-' + str(last_day)
    
    # Creating a tempdir to store the HDF
    tp = tempfile.mkdtemp()
    path_hdf = os.path.join(tp, 'HDF')
    os.makedirs(path_hdf, exist_ok = True)
        
    # create the connection with the https website
    down = pymodis.downmodis.downModis(path_hdf, password = pwd, user = usr, 
    url = "https://e4ftl01.cr.usgs.gov/", tiles = tile, path = path,
    product = product, today = last_date, enddate = first_date, debug = True)
    
    # initiate to the website download all the HDF for the month in path_hdf
    down.connect()
    down.downloadsAllDay()
#
    # Compute a raster with high quality pixel only from the HDF files.
    list_path_hdf = os.listdir(path_hdf)
    hdf_files = [f for f in list_path_hdf if ('.xml' not in f and 'hdf' in f)]
    list_lst_high = []
    list_non_missing = []
    for hdf_num in range(0, len(hdf_files)):
        hdf_ds = gdal.Open(os.path.join(path_hdf, hdf_files[hdf_num]), gdal.GA_ReadOnly)
        subdatasets = hdf_ds.GetSubDatasets()
    #
        #Reading the LST layer as an array
        lst_sds_name = subdatasets[lst_sds_id][0]
        lst_sds = gdal.Open(lst_sds_name, gdal.GA_ReadOnly)
        lst_array = lst_sds.ReadAsArray()
    #
        #Reading the QC layer as an array
        QC_sds_name = subdatasets[qc_sds_id][0]
        QC_sds = gdal.Open(QC_sds_name, gdal.GA_ReadOnly)
        QC_array = QC_sds.ReadAsArray()
    #
        # Extracting the georeferencing information (only needed once)
        transfo = lst_sds.GetGeoTransform()
        proj = lst_sds.GetProjection()
        #
        hdf_ds, lst_sds, QC_sds = None, None, None
        #
        LST_high = raster_high_qc(LST = lst_array, QC = QC_array)
        list_lst_high.append(LST_high[0])
        list_non_missing.append(LST_high[1]) 
    '''
    Averaging the dates to create a monthly average.
    In order to deal with missing values, the sum of the matrices for the month
    is divided by the matrix containing he number of valid observation 
    this month. This is done because numpy has trouble dealing with NaNs
    '''
    monthly_lst = sum(list_lst_high) / sum(list_non_missing)
    # replacing the nan (no observations at all) by 0
    monthly_lst[np.isnan(monthly_lst)] = 0 
    # switching from float to unsigned 16 bits integer to reduce the memory needed
    monthly_lst = monthly_lst.astype(np.dtype('u2'))
    name_geotiff = product +'.A' + year + month + '.' + tile + '.tif'
#
#
    # Save the monthly composite in a temporary folder to be used later on
    save_raster(value = monthly_lst, transfo = transfo, 
    projection = proj, save_path = os.path.join(temp_path_geotiff, name_geotiff), 
    in_na_value = 0)
    # Erasing the temporary hdf folder
    rmtree(path_hdf, ignore_errors = True)



'''
monthly_mosaic_lst
Creates a monhtly mosaic of several mosaic for a given month and a combination 
of parameters (product * year * month * day_night).
It computes a monthly temporary GeoTiff running df_month_hdf() for each tile, 
mosaics them, and save the resulting monthly mosaic as a GeoTiff.
The computation of individual monthly composite can be parallelized.
The temporary GeoTiff can be erased from the hard drive.
The mosaic is saved using the following convention 
'xxxx.Ayyyymm.day_night.tif' (xxxx: product name, yyyy: year,
mm: month', , day_night: time of the day)

Parameters:
product (str): Product to dowload, either 'MOD11A2.006' or 'MYD11A2.006'
day_night (str): Time of the day, either 'day' or 'night'
year (int): year for which the data must be processed
month (int): month for which the composite must be compiled
tilesH (list): A list of the horizontal tiles identifiers,
tilesV (list): A list of the vertical identifiers (e.g if one wants to 
    mosaic h00v02, h00v03, h00v03, tilesH = [0,0,0] and tilesV = [2,3,4].
save_path (str): path to save the mosaic. 
    Will be created if it does not exists.
usr (str): user name as entered for https://e4ftl01.cr.usgs.gov
pwd (str): password as entered for https://e4ftl01.cr.usgs.gov
ncores (int): If ncores > 1, multicore threading will be 
    used on a single machine.
'''
def monthly_mosaic_lst(product, month, year, day_night, tilesH,
    tilesV, save_path, usr, pwd, ncores):
    # Creating a tempdir to store the monthly geotiff per tile
    mosaic_name = (product +'.A' + str(year) + 
    str(month).zfill(2) + '.' + day_night)
    tp = tempfile.mkdtemp()
    temp_path_geotiff = os.path.join(tp, 'GEOTIFF', mosaic_name)
    os.makedirs(temp_path_geotiff, exist_ok = True)
    #
    # One monthly composite is produced for each tile (on one or several cores)
    Parallel(n_jobs = ncores, backend = 'threading')(
        delayed(dl_month_hdf)(product = product, day_night = day_night,
                tileH = tilesH[i], tileV = tilesV[i], year = year, 
                month = month, usr = usr, pwd = pwd,
                temp_path_geotiff = temp_path_geotiff) 
        for i in range(len(tilesV)))
#
    # Build a virtual raster mosaicing all the GeoTiff in temp_path_geotiff
    # VRT are saved in a temporary dir
    path_vrt = tempfile.mkdtemp()
    os.makedirs(path_vrt, exist_ok = True)
    VRT = gdal.BuildVRT(os.path.join(path_vrt, mosaic_name + ".vrt"),
    glob.glob(os.path.join(temp_path_geotiff, "*.tif")))

    # Save the mosaic as a GeoTiff
    save_raster(value = VRT.ReadAsArray(), transfo = VRT.GetGeoTransform(), 
    projection = VRT.GetProjection(), 
    save_path = os.path.join(save_path, mosaic_name + '.tif'), in_na_value = 0)
    VRT = None
    # Removing the temporary geotiff
    rmtree(temp_path_geotiff, ignore_errors = True)
    # Removing the VRT temporary directory
    rmtree(path_vrt, ignore_errors = True)


'''
mosaics
Wrap up function that creates multiple mosaics defined as combinations of product, year, month,
time of the day

Parameters:
parameters (list). Contains the required combinations of parameter.
    Each element is a list [product, year, month, day_night]
    Can be produced by mosaic_parameters()
tiles (list). Contains all the possible tiles identifiers.
    Each list is one dimensions [[list tilesH], [list tilesV]]
    Can be produced by existing_tiles()
ncores (int). If ncores > 1, multicore threading will be 
    used on a single machine.
save_path (str). A path pointing to the directory where the mosaic geotiff s
    hould be saved.
'''
def mosaics(save_path, tiles, parameters, usr, pwd, ncores):
    for unique_param in parameters:
        product_unique = unique_param[0]
        year_unique = unique_param[1]
        month_unique = unique_param[2]
        day_night_unique = unique_param[3]
#
        monthly_mosaic_lst(product = product_unique, month = month_unique, 
            year = year_unique, day_night = day_night_unique,
            tilesH = tiles[0], tilesV = tiles[1], 
            save_path = save_path, usr =usr, pwd = pwd, ncores = ncores)