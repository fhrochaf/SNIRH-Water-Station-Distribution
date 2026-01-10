import folium
import geopandas as gpd
from rasterio.enums import Resampling
import matplotlib as mpl
from matplotlib.colors import rgb2hex
import matplotlib.colors as mcolors
import numpy as np
import branca.colormap as cm
import pandas as pd

def display_Folium_map(geojson_display_dicts, zoom_start: int = 6):
    """
    Creates a Folium map from a GeoDataFrame, with hover and popup attributes.
    It will center on the first GeoDataFrame of the list

    geojson_display_dict = dictionary with geodataframes mapped to keyworded arguments to be passed to
    folium.GeoJson(...)

    zoom_start = Initial zoom level for the map.
    """

    print('Creating a Folium Map visualization.')

    # Define CRS projection to WGS 84 (Folium Standard)
    for i, dict in enumerate(geojson_display_dicts):
        if dict['data'].crs != 4326:
            dict['data'] = dict['data'].to_crs(epsg=4326)

    # Calculate map center
    center = [
        geojson_display_dicts[0]['data'].geometry.centroid.y.mean(),
        geojson_display_dicts[0]['data'].geometry.centroid.x.mean()
    ]

    # Initialize map
    fmap = folium.Map(
        location=center,
        zoom_start=zoom_start,
        tiles="OpenStreetMap"
    )

    for i in range(len(geojson_display_dicts)):
        #Check if the layer is of point type
        is_point_layer = geojson_display_dicts[i]['data'].geometry.iloc[0].geom_type == "Point"

        # Tooltip (hover)
        if 'attribute_map' in geojson_display_dicts[i]:
            tooltip = folium.GeoJsonTooltip(
                fields=list(geojson_display_dicts[i]['attribute_map'].keys()),
                aliases=[f"{label}:" for label in geojson_display_dicts[i]['attribute_map'].values()],
                localize=True
            )

            # Popup (click)
            popup = folium.GeoJsonPopup(
                fields=list(geojson_display_dicts[i]['attribute_map'].keys()),
                aliases=[f"{label}:" for label in geojson_display_dicts[i]['attribute_map'].values()],
                localize=True
            )
        else:
            tooltip = None
            popup = None

        # Define style function that allows: dict styles (callables), function styles (through lambda functions)
        if 'feature_settings' in geojson_display_dicts[i]:
            if callable(geojson_display_dicts[i]['feature_settings']):
                style = geojson_display_dicts[i]['feature_settings']
            else:
                style = lambda feature, i=i: geojson_display_dicts[i]['feature_settings']
        
        else: style = None

                  
        # Add geospatial data to map
        folium.GeoJson(
            geojson_display_dicts[i]['data'],
            name=geojson_display_dicts[i].get('name', f'Layer {i}'),  # Add this line
            zoom_on_click=True,
            tooltip=tooltip,
            marker=geojson_display_dicts[i].get('marker', None),
            popup=popup,
            style_function=None if is_point_layer else style,
            highlight_function=lambda feature: {
                "weight": 2,
                "color": "#ff7800"
                }
        ).add_to(fmap)

    folium.LayerControl().add_to(fmap)

    return fmap


def choose_subbasin(SUB_BASINS_PATH):
    """
    Ask for the user to choose for which one the data will be retrieved and analysis will be performed
    """
    gdf_subbasins = gpd.read_file(SUB_BASINS_PATH)
    valid_ids  = gdf_subbasins['DNS_NU_SUB'].astype(int)

    while True:
        print('Sub-basins: IDs\n')
        print(dict(zip(gdf_subbasins['DNS_NM'], valid_ids)))
        ID_SUBBASIN = input("\nType the code of the sub-basin (q to quit): ")

        if ID_SUBBASIN.lower() == "q":
            return None, None

        print('id subbasin', int(ID_SUBBASIN))
        print(type(int(ID_SUBBASIN)))

        try:
            if int(ID_SUBBASIN) in valid_ids.values:
                break
        except ValueError:
            pass
        print("Invalid input. Please try again.")

    return  gdf_subbasins[gdf_subbasins['DNS_NU_SUB'] == int(ID_SUBBASIN)], ID_SUBBASIN



def prepare_raster(raster, target_crs="EPSG:3857", resolution=120, method="nearest"):
    """Reproject and resample raster with method selection"""
    method_map = {
        "nearest": Resampling.nearest,
        "bilinear": Resampling.bilinear,
        "cubic": Resampling.cubic
    }

    try:
        raster =  raster.rio.reproject(
            target_crs, 
            resolution=resolution,
            resampling=method_map.get(method, Resampling.nearest)
        )
        print('Raster resampled successfully.')
    except Exception as e:
        print(f"Error resampling raster.")
        raise e
    
    return raster


def analyse_station_density_catchments(gdf_stations : gpd.GeoDataFrame,
                                       gdf_subbasin : gpd.GeoDataFrame,
                                       gdf_catchments : gpd.GeoDataFrame,
                                       gdf_streams : gpd.GeoDataFrame,
                                       min_strahler_order : int = 6):

    """
    Analyse the spatial density of monitoring stations within catchment areas and
    visualize the results together with hydrological features on an interactive map.

    This function:
    - Counts the number of monitoring stations located within each catchment polygon
    - Adds the station count as a new attribute to the catchments GeoDataFrame
    - Filters stream geometries by a minimum Strahler order
    - Creates a Folium map displaying stations, sub-basins, catchments, and streams
    with styling based on station presence and Strahler order

    Parameters
    ----------
    gdf_stations : geopandas.GeoDataFrame
        GeoDataFrame containing point geometries of monitoring stations.
        Must be in the same coordinate reference system as the other GeoDataFrames.

    gdf_subbasin : geopandas.GeoDataFrame
        GeoDataFrame containing polygon geometries of sub-basins to be displayed
        as contextual hydrological boundaries.

    gdf_catchments : geopandas.GeoDataFrame
        GeoDataFrame containing polygon geometries of catchment areas.
        A new column named 'count_st_cathcm' will be created to store the number
        of stations within each catchment.

    gdf_streams : geopandas.GeoDataFrame
        GeoDataFrame containing line geometries of river or stream networks.
        Must include a 'strord' column representing the Strahler stream order.

    min_strahler_order : int, optional
        Minimum Strahler order threshold used to filter which streams are displayed
        on the map. Default is 6.

    Returns
    -------
    gdf_catchments : geopandas.GeoDataFrame
        The input catchments GeoDataFrame with an additional column
        ('count_st_cathcm') representing the number of stations per catchment.

    map : folium.Map
        An interactive Folium map displaying stations, sub-basins, catchments,
        and filtered streams with dynamic styling.
    """
    

    # Count number of stations within each catchment area
    gdf_catchments["count_st_cathcm"] = (
        gpd.sjoin(gdf_stations,
                gdf_catchments, predicate="within")
        .groupby("index_right")
        .size()
        .reindex(gdf_catchments.index, fill_value=0)
    )


    ### Display results on a Folium MAP - stations and streams styled according to Strahler Order
    
    # Build a colormap for Streams, based on Strahler order
    min_order = gdf_streams["strord"].min()
    max_order = gdf_streams["strord"].max()
    strord_cmap = cm.linear.Blues_09.scale(min_order, max_order)
    strord_cmap.caption = "Strahler Order"

    # Style function for the streams, displayed according to Strahler order
    def stream_style_function(feature):
        order = feature["properties"]["strord"]

        return {
            "color": strord_cmap(order),
            "weight": order * 0.6,   # scale thickness
            "opacity": 0.9
        }


    # Build a colormap for catchments, based in the existence or absence of stations on them
    catch_cmap = cm.linear.RdYlGn_08.scale(0, 1)
    catch_cmap.caption = "Stations on catchments"

    # Style function for the catchments, displayed according to either having telemetering stations or not
    def catchment_style_function(feature):
        has_stations = bool(feature["properties"]["count_st_cathcm"])

        return {
            "fillColor": catch_cmap(has_stations), # Color the catchments according to the presence or absence of stations
            "color": "black",
            "weight": 1,
            "fillOpacity": 0.5,
        }
    
    # Filter which minimum stream Strahler Order will be displayed
    gdf_streams = gdf_streams[gdf_streams["strord"] >= min_strahler_order]

    geojson_display_dicts = [
        {
            'data' : gdf_stations,
            'attribute_map':
            {
                "codigoestacao": "Station ID",
                "Operadora_Sigla": "Responsable",
                "Municipio_Nome": "Municipality"
            },
            'feature_settings': {} #Point features, no style_function
        },
        {
            'data' : gdf_subbasin,
            'attribute_map' :
            {
                "DNS_DNB_CD": "Code of the basin",
                "DNS_NU_SUB": "Code of the sub-basin",
                "DNS_NM": "Name"
            },
            'feature_settings':
            {
                "color": "red",
                "weight": 3,
                "fill": False
            }, # Sub-basin feature style
        },
        {
            'data' : gdf_catchments,
            'attribute_map' :
            {
                "value": "#",
                "area_km2": "Area",
                "count_st_cathcm": "Station Count"
            },
            'feature_settings' : catchment_style_function, # Function to dinamically adjust the color of catchments 
        },
        {
            'data' : gdf_streams,
            'attribute_map' :
            {
                "strord": 'Strahler Order'
            },
            'feature_settings' : stream_style_function, # Function to dinamically adjust width and color of streams absed on order 
        }
    ]

    # Plots both the filtered stations with data, and the respective sub-basin
    map = display_Folium_map(geojson_display_dicts)

    return gdf_catchments, map


def analyse_station_upstream_network(river_net,
                                     gdf_subbasin,
                                     gdf_stations,
                                     gdf_catchments,
                                     gdf_streams
                                     ):
    """
    Analyse and visualize the upstream river network for each monitoring station.

    The function:
    - Extracts upstream river sub-networks for each station in a river network graph
    - Computes upstream stream lengths per Strahler order and total upstream length
    - Stores computed lengths as new attributes in the stations GeoDataFrame
    - Displays stations, catchments, sub-basins, and upstream networks on a Folium map

    Parameters
    ----------
    river_net : object
        River network object containing stream topology and station connectivity.

    gdf_subbasin : geopandas.GeoDataFrame
        GeoDataFrame containing sub-basin polygon geometries.

    gdf_stations : geopandas.GeoDataFrame
        GeoDataFrame containing station point geometries snapped to the river network.

    gdf_catchments : geopandas.GeoDataFrame
        GeoDataFrame containing catchment polygon geometries.

    gdf_streams : geopandas.GeoDataFrame
        GeoDataFrame containing stream geometries with Strahler order information.

    Returns
    -------
    gdf_stations : geopandas.GeoDataFrame
        Stations GeoDataFrame with added upstream length attributes.

    map : folium.Map
        Interactive map displaying the upstream networks and hydrological context.
    """


    def create_folium_colormap(base_color_hex, min_val, max_val):
        """
        Create a Folium colormap from light to base color.
        Used to create a different colormap for each station
        """
        base_rgb = mcolors.hex2color(base_color_hex)
        
        # Create color gradient
        colors = [
            rgb2hex(tuple(np.array(base_rgb) * 0.3 + np.array([1, 1, 1]) * 0.7)),
            rgb2hex(tuple(np.array(base_rgb) * 0.6 + np.array([1, 1, 1]) * 0.4)),
            base_color_hex,
            rgb2hex(tuple(np.array(base_rgb) * 0.8))
        ]
        
        # Create Folium colormap
        colormap = cm.LinearColormap(
            colors=colors,
            vmin=min_val,
            vmax=max_val,
            caption=f'Strahler Order'
        )
        
        return colormap

    # Check the stations belonging to the river network Graph
    station_ids = river_net.list_all_station_ids()

    geojson_display_dicts = [
        {
            'data' : river_net.streams_split,
            'name' : 'Streams',
            'attribute_map' : {},
            'feature_settings' : 
            {
                "color": "darkblue",
                "weight": 0.5,
                "alpha": 0.5
            }
        },
        {
            'data' : gdf_stations,
            'name' : 'Stations',
            'attribute_map' :
            {
                "codigoestacao": "Station ID",
                "Operadora_Sigla": "Responsable",
                "Municipio_Nome": "Municipality"
            },
            'feature_settings' : {}
        },
        {
            'data' : gdf_subbasin,
            'attribute_map' :
            {
                "DNS_DNB_CD": "Code of the basin",
                "DNS_NU_SUB": "Code of the sub-basin",
                "DNS_NM": "Name"
            },
            'feature_settings':
            {
                "color": "red",
                "weight": 3,
                "fill": False
            }, # Sub-basin feature style
        },
        {
            'data' : gdf_catchments,
            'name' : 'Sub catchments',
            'attribute_map' :
            {
                "value": "#",
                "area_km2": "Area"
            },
            'feature_settings' : 
            {
                "color": "black",
                "weight": 1,
                "fill": False
            }, # Function to dinamically adjust the color of catchments
        }
    ]

    # Create a base hue for each station sub-graph, based on Viridis colormap
    viridis = mpl.colormaps['viridis']
    station_colors_hex = {
        station_id: rgb2hex(viridis(i / len(station_ids)))
        for i, station_id in enumerate(station_ids)
    }

    min_order = gdf_streams["strord"].min()
    max_order = gdf_streams["strord"].max()

    station_colormaps = {
        station: create_folium_colormap(station_colors_hex[station], min_order, max_order)
        for station in station_ids
        }

    # Create empty geodataframe that will collect each influence zone
    gdf_station_upstreams = gpd.GeoDataFrame()

    # Create a display setting for each subgraph
    for station in station_ids:
        
        # Style function for the streams, displayed according to Strahler order
        def stream_style_fun_subgraph(feature):
            order = feature["properties"]["strahler"]
            station = feature["properties"]["dwnstream_stat"]
            strord_cmap = station_colormaps[station]

            return {
                "color": strord_cmap(order),
                "weight": order * 0.6,
                "opacity": 0.6
            }
            
        gdf_subgraph = river_net.get_upstream_geodataframe(station,
                                                        include_station_node=True,
                                                        include_nodes=False)
        gdf_subgraph['dwnstream_stat'] = station

        # Add display settings of the subgraph to the display dict
        geojson_display_dicts.append({
            'data' : gdf_subgraph,
            'name' : f'Upstream network of station {station}',
            'attribute_map' :
            {
                "strahler": "Strahler Order",
                "dwnstream_stat" : "Downstream station"
            },
            'feature_settings' : stream_style_fun_subgraph # Function to dinamically adjust width and color of streams absed on order
        })

        # Add gdf_subgraph to the geodataframe to be returned
        gdf_station_upstreams = pd.concat([gdf_subgraph, gdf_station_upstreams])


    # Plots both the filtered stations with data, and the respective sub-basin
    map = display_Folium_map(geojson_display_dicts)

    return gdf_station_upstreams, map


def analyse_station_influence_zone(river_net,
                                   gdf_subbasin,
                                   gdf_stations,
                                   gdf_catchments,
                                   gdf_streams
                                   ):
    """
        Analyse and visualize the influence zones of monitoring stations within a river network.

        The function:
        - Computes upstream influence zones for each station in the river network
        - Visualizes stations, sub-basins, catchments, streams, and influence zones
        on an interactive Folium map, styled by Strahler order

        Parameters
        ----------
        river_net : object
            River network object containing stream topology and station connectivity,
            with methods to compute and retrieve station influence zones.

        gdf_subbasin : geopandas.GeoDataFrame
            GeoDataFrame containing sub-basin polygon geometries.

        gdf_stations : geopandas.GeoDataFrame
            GeoDataFrame containing station point geometries snapped to the river network.

        gdf_catchments : geopandas.GeoDataFrame
            GeoDataFrame containing catchment polygon geometries.

        gdf_streams : geopandas.GeoDataFrame
            GeoDataFrame containing stream geometries with Strahler order information.

        Returns
        -------
        map : folium.Map
            Interactive map displaying station influence zones and hydrological context.

        gdf_influence_zones : geopandas.GeoDataFrame
            GeoDataFrame containing stream geometries upstream up to the maximum influence zone of the stations.
        """

    # Check the stations belonging to the river network Graph
    station_ids = river_net.list_all_station_ids()

    # Compute influence zones (e.g., 100 km upstream from each station)
    river_net.compute_station_influence_zones()

    geojson_display_dicts = [
        {
            'data' : river_net.streams_split,
            'name' : 'Streams',
            'attribute_map' : {},
            'feature_settings' : 
            {
                "color": "darkblue",
                "weight": 0.5,
                "alpha": 0.5
            }
        },
        {
            'data' : gdf_stations,
            'name' : 'Stations',
            'attribute_map' :
            {
                "codigoestacao": "Station ID",
                "Operadora_Sigla": "Responsable",
                "Municipio_Nome": "Municipality"
            },
            'feature_settings' : {}
        },
        {
            'data' : gdf_subbasin,
            'attribute_map' :
            {
                "DNS_DNB_CD": "Code of the basin",
                "DNS_NU_SUB": "Code of the sub-basin",
                "DNS_NM": "Name"
            },
            'feature_settings':
            {
                "color": "red",
                "weight": 3,
                "fill": False
            }, # Sub-basin feature style
        },
        {
            'data' : gdf_catchments,
            'name' : 'Sub catchments',
            'attribute_map' :
            {
                "value": "#",
                "area_km2": "Area"
            },
            'feature_settings' : 
            {
                "color": "black",
                "weight": 1,
                "fill": False
            }, # Function to dinamically adjust the color of catchments
        }
    ]

    # Build a colormap for Streams, based on Strahler order
    min_order = gdf_streams["strord"].min()
    max_order = gdf_streams["strord"].max()
    strord_cmap = cm.linear.RdPu_09.scale(min_order, max_order)
    strord_cmap.caption = "Strahler Order"

    # Create empty geodataframe that will collect each influence zone
    gdf_influence_zones = gpd.GeoDataFrame()

    # Style function for the streams, displayed according to Strahler order
    def stream_style_function(feature):
        order = feature["properties"]["strahler"]

        return {
            "color": strord_cmap(order),
            "weight": order * 0.6,   # scale thickness
            "opacity": 0.9
        }


    # Create a display setting for each subgraph
    for station in station_ids:
        
    # Get influence zone for specific station
        station_zone_gdf = river_net.get_station_influence_geodataframe(station)

        # Add display settings of the subgraph to the display dict
        geojson_display_dicts.append({
            'data' : station_zone_gdf,
            'name' : f'{station} station recommended influence zone.',
            'attribute_map' :
            {
                "strahler": "Strahler Order",
                "closest_station": "Closest Station"
            },
            'feature_settings' : stream_style_function # Function to dinamically adjust width and color of streams absed on order
        })

        gdf_influence_zones = pd.concat([gdf_influence_zones, station_zone_gdf])

    # Plots both the filtered stations with data, and the respective sub-basin
    map = display_Folium_map(geojson_display_dicts)

    return gdf_influence_zones, map