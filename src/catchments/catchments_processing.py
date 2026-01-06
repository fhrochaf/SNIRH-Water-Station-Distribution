import os

import geopandas as gpd
import numpy as np
from rasterio.features import rasterize
from rasterio import features
import pyflwdir
from shapely.ops import nearest_points
import pandas as pd
import rioxarray
import geopandas as gpd
import pandas as pd
from shapely.ops import linemerge, unary_union
import rasterio
from rasterio.io import MemoryFile
import requests
import branca.colormap as cm
from utils import *

def get_dem(bounds,
            asset_url : str = None,
            key_opentopo : str = None,
            save_dem_raster : bool = False,
            output_dir : str = None,
            file_name : str = 'dem'):
    
    """
    Get DEM data a final end point from a WMS Service.
    Standard = from Global Asset (OpenLandMap STAC), 30m resolution, cropped by the bounding box of your interest

    bounds = [xmin, ymin, xmax, ymax]
    asset_url = end point URL for a raster in a webMapService
    resample_spatial_res = desired spatial resolution, in meters
    """

    if output_dir is None and os.path.isdir(os.path.join(os.getcwd(), 'outputs')):
        output_dir = os.path.join(os.getcwd(), 'outputs')
    elif output_dir is None:
        os.mkdir(os.path.join(os.getcwd(), 'outputs'))
        output_dir = os.path.join(os.getcwd(), 'outputs')

    # # Asset: global DSM 30m resolution data (OpenLandMap STAC / Copernius DEM digital surface model / dsm_glo30_20110101_20151231)
    #     asset_url = "https://s3.openlandmap.org/arco/dsm_glo30_m_30m_s_20110101_20151231_go.epsg.4326_v20211004.tif"
 
    if not asset_url:
        asset_url = "https://s3.openlandmap.org/arco/dsm_glo30_m_30m_s_20110101_20151231_go.epsg.4326_v20211004.tif"


    # First try to Get DEM from Open Land Map
    print(f'Attempting to load DEM from {asset_url}.')
    try:
        data = rioxarray.open_rasterio(asset_url, masked=True)
        dem_data = data.rio.clip_box(*bounds)

        crs = dem_data.rio.crs
        min = dem_data.min()
        max = dem_data.max()

        if save_dem_raster:    
            print("Writing DEM file.")
            try:

                print(f"CRS: {crs}.")
                print(f"Data range: {float(min):.2f}m to {float(max):.2f}m.")

                file_path = os.path.join(output_dir, f'{file_name}.tif')

                dem_data.rio.to_raster(file_path)
                return dem_data
                    
            except Exception as e:
                print(f"Failed to write file: {e}")

    except Exception as e:
        print(f"Failed to get DEM on Open Land Map {asset_url}: {e}")

        if key_opentopo:
            # If failed, the DEM will be requested to Open Topography
            print("Attempting connection on Open Topography.")
            try:
                asset_url = 'https://portal.opentopography.org/API/globaldem'

                west, south, east, north = bounds.astype(str)
                params = {
                    "accept": "*/*",
                    "demtype": "SRTMGL3",
                    "south": south,
                    "north": north,
                    "west": west,
                    "east": east,
                    "outputFormat": "GTiff",
                    "API_Key": key_opentopo
                    }
            
                response = requests.get(asset_url, params=params)
                asset_url = response.url

                # data = rioxarray.open_rasterio(asset_url, masked=True)
                # dem_data = data.rio.clip_box(*bounds)
                
                print(f'Attempting to load DEM from {asset_url}.')
                with MemoryFile(response.content) as memfile:
                    with memfile.open() as src:
                        dem_data = src.read()
                        dem_meta = src.meta.copy()
                        crs = src.crs
                        nodata = src.nodata
                        transform = src.transform
                min = dem_data.min()
                max = dem_data.max()
                
                dem_meta.update({
                "transform": transform,
                "driver":"GTiff",
                "height": dem_data.shape[1],
                "width": dem_data.shape[2],
                "nodata": nodata})

            except Exception as e:
                print(f'Failed to get DEM on Open Topography: {e}')
                return None
            

    if save_dem_raster:    
        print("Writing DEM file.")
        try:

            print(f"CRS: {crs}.")
            print(f"Data range: {float(min):.2f}m to {float(max):.2f}m.")

            file_path = os.path.join(output_dir, f'{file_name}.tif')

            with rasterio.open(file_path, "w", **dem_meta) as dst:
                dst.write(dem_data)
                
        except Exception as e:
            print(f"Failed to write file: {e}")

    return dem_data


# convenience method for vectorizing a raster
def vectorize(data, nodata, transform, crs, name="value"):
    """Vectorize feature using rasterio features vectorize
    """
    feats_gen = features.shapes(
        data,
        mask=data != nodata,
        transform=transform,
        connectivity=8,
    )
    feats = [
        {"geometry": geom, "properties": {name: val}} for geom, val in list(feats_gen)
    ]

    # parse to geopandas for plotting / writing to file
    gdf = gpd.GeoDataFrame.from_features(feats, crs=crs)
    gdf[name] = gdf[name].astype(data.dtype)
    return gdf



def pyflwdir_subbasins_minarea(
        dem,
        gdf : gpd.GeoDataFrame,
        min_sto : int = 4,
        save_data : bool = False,
        output_dir: str = None,
        file_name1 : str = 'catchments',
        file_name2 : str = 'streams',
        file_name3 : str = 'flowdir'):
    """
    Delineate Pfafstetter sub-basins and return a GeoDataFrame from pyflwdir library.
    dem = rioxarray loaded dem
    min_area = minimum area for the catchments, km²
    gdf = geopandas dataframe of the larger sub_basin
    depth = depth of sub_basin delineating scale (the higher, the smaller and more numerous are the basins)
    min_sto = Minimum Strahler Order recognized as river, by the default 4.
    """

    if output_dir is None and os.path.isdir(os.path.join(os.getcwd(), 'outputs')):
        output_dir = os.path.join(os.getcwd(), 'outputs')
    elif output_dir is None:
        os.mkdir(os.path.join(os.getcwd(), 'outputs'))
        output_dir = os.path.join(os.getcwd(), 'outputs')


    # Compute total area of the GeoDataFrame (km²)
    total_area_km2 = gdf.to_crs(3857).geometry.area.sum() / 1e6

    while True:
        try:
            min_area = input(
                f"Choose the minimum area to delineate catchment areas "
                f"(in km², max {int(total_area_km2)}). Type q to cancel and quit: "
            )

            # Check empty input
            if min_area.lower() == 'q':
                return None, None
            
            # Check empty input
            if not min_area.strip():
                raise ValueError("Input cannot be empty.")

            min_area = int(min_area)

            # Check positivity
            if min_area <= 0:
                raise ValueError("Area must be a positive integer.")

            # Check against basin size
            if min_area > total_area_km2:
                raise ValueError(
                    f"Area cannot exceed basin area ({total_area_km2:.2f} km²)."
                )

            break  # ✅ valid input

        except ValueError as e:
            print(f"Invalid input: {e}")


    # Get metadata from raster
    crs = dem.rio.crs # This should be your target CRS
    width = dem.rio.width
    height = dem.rio.height
    count = dem.rio.count

    print(f"Raster CRS: {crs}")
    print(f"Input GeoDataFrame CRS: {gdf.crs}")
    
    try:
        # Reproject gdf to match the raster CRS
        if gdf.crs != crs:
            print(f"Reprojecting GeoDataFrame from {gdf.crs} to {crs}")
            gdf = gdf.to_crs(crs=crs)
        else:
            print("GeoDataFrame already in correct CRS")
    except Exception as e:
        print(f'Failed to reproject GeoDataFrame.')

    # Nodata value
    nodata = dem.rio.nodata
    # Get affine transform array
    transform = dem.rio.transform()
    # Squeeze raster if it has more than 2 dimensions (rioxarray tipically (1,x,y))
    if len(dem.shape) > 2:
        dem = dem.squeeze()

    # Get the values / transform into a numpy array
    dem = dem.values

    # Rasterize basin
    mask = rasterize(
        [(geom, 1) for geom in gdf.geometry],
        out_shape=dem.shape,
        transform=transform,
        fill=0,
        dtype="uint8"
    )

    # Mask DEM outside basin geometry
    dem_masked = np.where(mask == 1, dem, nodata)

    # Create flow direction array
    print('Calculating flow directions for the elevation data.')
    flwdir = pyflwdir.from_dem(
        dem_masked,
        transform=transform,
        nodata=nodata,
        latlon=False  # IMPORTANT for Pfafstetter
    )

    # Delineate basins according to minimum area given
    pfaf, _ = flwdir.subbasins_area(min_area)

    # Vectorize features

    gdf_catchments = vectorize(pfaf.astype(np.int32), 0,transform, crs, name="value")
    gdf_catchments['area_km2'] = gdf_catchments.area/(1e6)
    # Remove micro basins, below the specified minimum area
    gdf_catchments = gdf_catchments[gdf_catchments['area_km2'] >= 2000]
    print(f'{len(gdf_catchments)} catchments computed.')

    # Vectorize streams and create a gdf for it
    print('Vectorizing streams.')
    feats = flwdir.streams(min_sto=min_sto)
    gdf_streams = gpd.GeoDataFrame.from_features(feats, crs=crs)
    # Select only streams inside our final catchments
    gdf_streams = gdf_streams.clip(gdf_catchments)
    # Select only line/multiline geometries
    gdf_streams = gdf_streams[gdf_streams.geometry.geom_type.isin(['LineString', 'MultiLineString'])]

    # Write dataset to file
    if save_data:
        try:
            print(f"Saving catchments of basin.")
            gdf_catchments.to_file(os.path.join(output_dir, f'{file_name1}_{min_area}.geojson'), driver='GeoJSON')

            print(f"Saving streams of basin.")
            gdf_streams.to_file(os.path.join(output_dir, f'{file_name2}.geojson'), driver='GeoJSON')

            print(f"Saving flow direction raster.")
            # With the to_array() method, pyflwdir can return a flow direction numpy array from the FlwDirRaster object in any supported convention.
            d8_data = flwdir.to_array(ftype="d8")

            # update data type and nodata value properties which are different compared to the input elevation grid and write to geotif
            prof = {
                'driver': 'GTiff',
                'dtype': d8_data.dtype,
                'count': count,
                'nodata': 247,
                'width': width,
                'height': height,
                'crs': crs,
                'transform': transform
            } #Profile dictionary to save the Flow Direction raster with proper metadata
            file_path = os.path.join(output_dir, f'{file_name3}.tif')
            with rasterio.open(file_path, "w", **prof) as src:
                src.write(d8_data, 1)
            
        except Exception as e:
            print(f'Files could not be saved: {e}')

    return gdf_catchments, gdf_streams


def snap_stations(
    gdf_stations : gpd.GeoDataFrame,
    gdf_to_snap : gpd.GeoDataFrame):
    """ Snaps stations to the geometries of target GeoDataFrame

    Creates a 'snapped' flag column to show that the point had its original position adjusted.
    """
    # Enforce an equal CRS for the GeoDataFrames:
    target_crs = 3857

    points = gdf_stations.to_crs(target_crs)
    targets = gdf_to_snap.to_crs(target_crs)

    # Check stations that are outside of the basin
    targets_union = targets.unary_union
    points_inside = points[points.within(targets_union)]
    points_inside['snapped'] = False
    points_outside = points[~points.within(targets_union)]


    if len(points_outside) > 0:
        # Find neares geometry for the point outside
        nearest_geoms = []
        for p in points_outside.geometry:
            nearest_geoms.append(nearest_points(p, targets_union)[1]) #[1] to return the closest vertex of the geometry 2 (poly_union), which p will be snapped to
        points_outside = points_outside.copy()
        points_outside['geometry'] = nearest_geoms
        points_outside['snapped'] = True

    # Merge points outside and inside
    points_final = pd.concat([points_inside, points_outside], ignore_index=True)
    points_final = gpd.GeoDataFrame(points_final, crs=target_crs)

    print(f'Stations have been snapped to the target geometry: {len(points_outside)}.')

    return points_final



def connect_streams(gdf,
                    attribute_column='strord',
                    tolerance : float = 1e-4,
                    save_data : bool = False,
                    output_dir: str = None,
                    file_name : str = 'streams_conn'):
    """
    Connect streams that are connected and have the samer Strord (or other provided attribute) value,
    based on an index that informs the connection of the segments ('idx_ds') in the case of streams created
    through pyflwdir_subbasins_minarea()

    gdf_streams = GeoDataFrame with streams
    segm_attr = The index that informs the connected streams
    tolarance = The maximum distance that features may be merged together
    agg_column = The attribute which each connected stream segment will be aggregated on
    """

    # Check directories
    if output_dir is None and os.path.isdir(os.path.join(os.getcwd(), 'outputs')):
        output_dir = os.path.join(os.getcwd(), 'outputs')
    elif output_dir is None:
        os.mkdir(os.path.join(os.getcwd(), 'outputs'))
        output_dir = os.path.join(os.getcwd(), 'outputs')
        
    # Make a copy
    gdf = gdf.copy()
    gdf['_orig_index'] = range(len(gdf))
    
    # Build spatial index
    spatial_index = gdf.sindex
    
    # Dictionary to track merges
    merge_groups = {}
    
    print(f"Computing connected geometries.")
    # Find touching segments with same attribute
    for idx, row in gdf.iterrows():
        attr_value = row[attribute_column]
        geom = row.geometry
        
        # Get potential intersections using spatial index
        possible_matches_index = list(spatial_index.intersection(geom.bounds))
        possible_matches = gdf.iloc[possible_matches_index]
        
        for match_idx, match_row in possible_matches.iterrows():
            if match_idx <= idx:
                continue  # Avoid duplicate comparisons
            
            if match_row[attribute_column] != attr_value:
                continue  # Different attribute value
            
            if geom.touches(match_row.geometry) or geom.intersects(match_row.geometry):
                if idx not in merge_groups:
                    merge_groups[idx] = set([idx])
                    merge_groups[idx].add(match_idx)
    
    # Merge segments
    merged_segments = []
    processed = set()
    
    for group_idx, segment_indices in merge_groups.items():
        if group_idx in processed:
            continue
        
        # Get all segments in this merge group
        segments_to_merge = list(segment_indices)
        
        # Find all connected segments recursively
        to_process = segments_to_merge.copy()
        while to_process:
            current = to_process.pop()
            if current in merge_groups:
                for connected_idx in merge_groups[current]:
                    if connected_idx not in segments_to_merge:
                        segments_to_merge.append(connected_idx)
                        to_process.append(connected_idx)
        
        # Mark as processed
        processed.update(segments_to_merge)
        
        # Get geometries to merge
        geoms_to_merge = [gdf.iloc[idx].geometry for idx in segments_to_merge]
        
        # Merge geometries
        try:
            merged_geom = linemerge(unary_union(geoms_to_merge))
            
            # Get attributes from first segment
            first_row = gdf.iloc[segments_to_merge[0]].copy()
            first_row.geometry = merged_geom
            first_row['_merged_from'] = segments_to_merge
            
            merged_segments.append(first_row)
        except:
            # If merge fails, keep original segments
            for idx in segments_to_merge:
                merged_segments.append(gdf.iloc[idx])
    
    # Add unmerged segments
    for idx in range(len(gdf)):
        if idx not in processed:
            merged_segments.append(gdf.iloc[idx])
    
    # Create result GeoDataFrame
    result_gdf = gpd.GeoDataFrame(merged_segments, crs=gdf.crs)
    
    # Clean up
    result_gdf = result_gdf.drop(columns=['_orig_index', '_merged_from'], errors='ignore')
   # Filter out all Point geometries
    result_gdf = result_gdf[~result_gdf.geometry.type.isin(['Point', 'MultiPoint'])]

    # Write dataset to file
    if save_data:
        try:
            print(f"Saving connected geometries.")
            result_gdf.to_file(os.path.join(output_dir, f'{file_name}.geojson'), driver='GeoJSON')
        except Exception as e:
            print(f'File could not be saved: {e}')


    return result_gdf


def analyse_station_density_catchments(gdf_stations : gpd.GeoDataFrame,
                                       gdf_subbasin : gpd.GeoDataFrame,
                                       gdf_catchments : gpd.GeoDataFrame,
                                       gdf_streams : gpd.GeoDataFrame,
                                       min_strahler_order : int = 6): 
    

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