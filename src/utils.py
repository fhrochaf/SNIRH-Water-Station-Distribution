import folium
import geopandas as gpd
from rasterio.enums import Resampling


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