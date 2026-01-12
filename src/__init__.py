import os
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns

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