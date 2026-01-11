import os
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv

load_dotenv()

# User must provide the credentials for HidroWeb access in the config.env file.
ID_SNIRH = os.getenv('ID_SNIRH') 
PWD_SNIRH = os.getenv('PWD_SNIRH')

# User must provide the credentials for OpenTopography access in the config.env file
KEY_OPEN_TOPOGRAPHY = os.getenv('KEY_OPEN_TOPOGRAPHY') 
# Alternativelly, DEM's within the bounding-box of each brazilian sub-basin are available in
# https://drive.google.com/drive/folders/1c_NcoKtEFZkNS9EsjcDbT4IeXCWQCbMT?usp=drive_link

# Define paths
ROOT_DIR = os.path.dirname(__file__)
os.chdir(ROOT_DIR)
OUTPUT_DIR = os.path.join(ROOT_DIR, 'outputs')
SUB_BASINS_PATH = os.path.join(ROOT_DIR, 'data/hidrosubbasins.geojson')

# Import all modules
from . import utils
from .hidroweb import hidrowebANA
from .catchments import catchments_processing, river_network

# Export everything
from .utils import *
from .hidroweb.hidrowebANA import *
from .catchments.catchments_processing import *
from .catchments.river_network import *

# Helper function to load subbasins
def load_subbasins():
    """Load Brazilian sub-basins GeoDataFrame"""
    return gpd.read_file(SUB_BASINS_PATH)