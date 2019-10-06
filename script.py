#########################################################
#       Downloading and mosaicing MODIS products        #
#########################################################
'''
Tested on I7-5600QM @ 2.60GHz / Win/, 12 GB RAM with Python 3.6.6 in an Anaconda Environment

PREREQUISITES:
The following libraries must be installed in the Python 
environnment used (available through pip):
* tempfile
* pymodis
* osgeo
* os
* calendar
* numpy as np
* shutil
* csv
* urllib.request
* io
* datetime
* joblib
* glob

The final user should have a NASA Earth data login and password 
'''

'''
ALGORITHM:
The functions are commented in details (in the functions.py file), 
but the algorithm works with four main functions:
* raster_high_qc() takes a LST and a QC numpy array as input, and return an 
        array with only the highest quality pixels
* dl_monthly_hdf() download all the HDFs for a tile for a unique combination 
    of parameters (product * year * month * day_or_night), passes them to 
    raster_high_qc() to keep only the highest quality pixels, averages them
    to create a monthly composite of the tile and saves them as temporary GeoTiffs.
* monthly_mosaic_lst() creates a monthly mosaic of several tiles for a unique 
    combination of parameters (product * year * month * day_or_night) and 
    saves it as a GeoTiff
* mosaics() creates monthly mosaics of several tiles for several combinations 
    of parameters and saves them as GeoTiff

NOTES: 
The script can be use in multicore by building a tile monthly composite in each core.
Only one tile is computed on each core to avoid using too much RAM (and too large overhead).
HDF and temporary GeoTiff are removed from the hard drive after they have been used to 
avoid excessive HD usage.
Moreover, unique combinations of parameters and tiles identifiers can be generated 
automatically using functions described in the function.py file, meaning that this 
function can be run on a cluster using a bash script. The only variables to pass 
to each node are (1) a path to save the final raster,  (2) an integer between 1 and 
the number of unique combination of parameters (here 832, indicating to each core 
which mosaic could be compiled). Any other variables (user name, password...) can be 
hardcoded in the python script.
The script sometimes send the following warning message:
'RuntimeError: Failed to write .vrt file in FlushCache()'. Though it does not impair the
computation, changing the working directory night be a quick fix around the issue.
'''

#########################################################
#                    Running script_test                #
#########################################################
# Importing the functions
import sys
#working_dir = 'path/to/directory'
working_dir = 'C:/Users/anmi/Downloads/script_test'
sys.path.append(working_dir)

import functions as mod

# Creating lists of all existing tiles
tiles = mod.existing_tiles()

# Creating a list of all possbile parameters combination
parameters = mod.mosaics_parameters()

# Provinding login informations
usr = 'carwoo'
pwd = '_Zbca3u4u'
mosaic_path = 'path/to/save_directory'

mod.mosaics(save_path = mosaic_path, tiles = tiles,
    parameters = parameters, usr = usr, pwd = pwd, ncores = 1)

t = time.time()
mod.mosaics(save_path = 'C:/Mosaic', tiles = [[22,22,22], [3,4,5]],
    parameters = parameters[0:2], usr = usr, pwd = pwd, ncores = 3)
print(time.time() - t)
