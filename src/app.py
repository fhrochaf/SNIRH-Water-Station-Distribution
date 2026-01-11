import os
import matplotlib.pyplot as plt
import seaborn as sns

def initialize(target_globals):
    ##########################
    #  Define Global Variables and Paths
    ##########################

    # Functions implemented in hidrowebANA.py: get_access_token(), get_subbasins(),
    # get_telemstat_subbasin(), get_telemstat_flowseries(), get_telemstat_flowseries_all(), filter_and_save_successfull_flowseries(),
    # display_folium_map(), choose_subbasin()

    # Functions implemented in catchments_procesing.py: get_dem(), vectorize(), pyflwdir_subbasins_minarea(),
    # snap_stations_to_subbasin() and connect_streams()

    ##########################
    #  Define Global Variables and Paths
    ##########################

    # Define the working directory as the same of this python file
    

    target_globals['ROOT_DIR'] = os.path.join(os.path.dirname(__file__))
    os.chdir(target_globals['ROOT_DIR'])
    print("Scripting running on current directory: ", os.getcwd())

    target_globals['OUTPUT_DIR'] = os.path.join(target_globals['ROOT_DIR'], 'outputs')

    # Define path of local needed data
    target_globals['SUB_BASINS_PATH'] = os.path.join(target_globals['ROOT_DIR'], 'data/hidrosubbasins.geojson')

    # Import functions, libraries and variables in the Utils module
    import utils
        # Add all public attributes to global namespace
    for name in dir(utils):
        if not name.startswith('_'):
            target_globals[name] = getattr(utils, name)


    # Import functions, libraries and variables in the Hidroweb module
    import hidroweb.hidrowebANA
        # Add all public attributes to global namespace
    for name in dir(hidroweb.hidrowebANA):
        if not name.startswith('_'):
            target_globals[name] = getattr(hidroweb.hidrowebANA, name)


    # Import functions, libraries and variables in the catchments_processing module
    import catchments.catchments_processing
        # Add all public attributes to global namespace
    for name in dir(catchments.catchments_processing):
        if not name.startswith('_'):
            target_globals[name] = getattr(catchments.catchments_processing, name)


    # Import functions, libraries and variables in the catchments_processing module
    import catchments.river_network
        # Add all public attributes to global namespace
    for name in dir(catchments.river_network):
        if not name.startswith('_'):
            target_globals[name] = getattr(catchments.river_network, name)



