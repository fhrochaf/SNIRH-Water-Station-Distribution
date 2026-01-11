import os

import geopandas as gpd
import numpy as np
from rasterio.features import rasterize
from rasterio import features
import pyflwdir
from shapely.ops import nearest_points
import pandas as pd
import geopandas as gpd
import pandas as pd
from shapely.ops import linemerge, unary_union
import rasterio
from rasterio.io import MemoryFile
import rioxarray
import requests

import rioxarray


from rasterio.enums import Resampling

def get_dem(bounds, key_opentopo=None, 
            save_dem_raster=False, output_dir=None, file_name='dem'):
    """
    Get DEM data from a web service, with fallback to OpenTopography.

    bounds: [xmin, ymin, xmax, ymax] bounding box
    key_opentopo: API key for OpenTopography
    save_dem_raster: Whether to save the DEM as a GeoTIFF
    output_dir: Output directory (defaults to ./outputs)
    file_name: Output filename without extension
        
    Returns:
        tuple: (dem_data, dem_meta) or None if failed
    """
    
    # Setup output directory
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'outputs')
        os.makedirs(output_dir, exist_ok=True)
    
    # Attempt to retrieve data from OpenTopography 
    dem_data = _load_from_opentopo(bounds, key_opentopo)
    
    if dem_data is None:
        print("Failed to retrieve DEM from all sources.")
        return None
    
    # Print metadata
    crs = dem_data.rio.crs
    data_min, data_max = dem_data.min(), dem_data.max()
    print(f"CRS: {crs}")
    print(f"Data range: {data_min:.2f}m to {data_max:.2f}m")
    
    dem_data = dem_data.squeeze() #Squeeze the raster into shape (x,y)
    # Save if requested
    if save_dem_raster:
        file_path = os.path.join(output_dir, f'{file_name}.tif')
        """Save raster data to file."""
        print(f"Saving DEM to {file_path}")
        try:
            dem_data.rio.to_raster(file_path)
        except Exception as e:
            print(f"Failed to save file: {e}")
 
    return dem_data


def _load_from_opentopo(bounds, api_key, target_crs="EPSG:3857", resolution=120, method="nearest"):
    """Load DEM from OpenTopography API and return as rioxarray with target CRS and resolution."""
    print('Attempting to load DEM from OpenTopography.')
    try:
        west, south, east, north = map(str, bounds)
        params = {
            "demtype": "SRTMGL3",
            "south": south,
            "north": north,
            "west": west,
            "east": east,
            "outputFormat": "GTiff",
            "API_Key": api_key
        }
        
        method_map = {
            "nearest": Resampling.nearest,
            "bilinear": Resampling.bilinear,
            "cubic": Resampling.cubic,
            "average": Resampling.average
        }

        url = 'https://portal.opentopography.org/API/globaldem'
        response = requests.get(url, params=params)
        response.raise_for_status()
        
        with MemoryFile(response.content) as memfile:
            with memfile.open() as src:
                # Read directly as rioxarray
                dem_raster = rioxarray.open_rasterio(src, masked=True)
                # Load into memory to avoid issues when MemoryFile closes
                dem_raster = dem_raster.load()

                # Reproject and resample to target CRS and resolution
                print('Resampling and reprojecting raster.')
                dem_raster = dem_raster.rio.reproject(
                    target_crs,
                    resolution=resolution,
                    resampling=method_map.get(method, Resampling.bilinear)
                )
                
        return dem_raster
    
    except Exception as e:
        print(f'Failed to load from OpenTopography: {e}')
        return None


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
    output_dir: Output directory (defaults to ./outputs)
    file_name: Output filename without extension
    """

    # Setup output directory
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'outputs')
        os.makedirs(output_dir, exist_ok=True)


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

            break  # valid input

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


def snap_stations_polyg(
    gdf_stations: gpd.GeoDataFrame,
    gdf_to_snap: gpd.GeoDataFrame,
    max_distance: float = None
) -> gpd.GeoDataFrame:
    """Snaps stations to the nearest geometry in target GeoDataFrame.
    
    Parameters
    ----------
    gdf_stations : gpd.GeoDataFrame
        GeoDataFrame containing station points to snap
    gdf_to_snap : gpd.GeoDataFrame
        GeoDataFrame containing target geometries to snap to
    max_distance : float, optional
        Maximum distance (in target CRS units) for snapping. 
        Points beyond this distance won't be snapped. If None, all points are snapped.
    
    Returns
    -------
    gpd.GeoDataFrame
        Stations with snapped geometry and additional columns:
        - 'snapped': boolean flag indicating if point was moved
        - 'snap_distance': distance the point was moved (in target CRS units)
    """
    # Enforce an equal CRS for the GeoDataFrames
    target_crs = 3857
    
    points = gdf_stations.to_crs(target_crs).copy()
    targets = gdf_to_snap.to_crs(target_crs)
    
    # Create union of all target geometries
    targets_union = targets.unary_union
    
    # Find nearest point on target geometries for each station
    snapped_geoms = []
    snap_distances = []
    
    for point in points.geometry:
        nearest_point = nearest_points(point, targets_union)[1]
        snapped_geoms.append(nearest_point)
        snap_distances.append(point.distance(nearest_point))
    
    # Create output GeoDataFrame
    points['snap_distance'] = snap_distances
    points['original_geometry'] = points.geometry
    points['geometry'] = snapped_geoms
    
    # Apply max_distance threshold if specified
    if max_distance is not None:
        points.loc[points['snap_distance'] > max_distance, 'geometry'] = points.loc[
            points['snap_distance'] > max_distance, 'original_geometry'
        ]
        points['snapped'] = points['snap_distance'] <= max_distance
        num_snapped = points['snapped'].sum()
    else:
        points['snapped'] = points['snap_distance'] > 0
        num_snapped = points['snapped'].sum()
    
    # Drop the temporary original_geometry column
    points = points.drop(columns=['original_geometry'])
    
    print(f'Stations snapped: {num_snapped} out of {len(points)}')
    
    return points


def snap_stations_streams(
    gdf_stations: gpd.GeoDataFrame,
    gdf_to_snap: gpd.GeoDataFrame,
    tolerance: float = 150  # snapping tolerance in meters
) -> gpd.GeoDataFrame:
    """Snaps stations to the nearest geometry in target GeoDataFrame using Shapely's snap.
    
    Parameters
    ----------
    gdf_stations : gpd.GeoDataFrame
        GeoDataFrame containing station points to snap
    gdf_to_snap : gpd.GeoDataFrame
        GeoDataFrame containing target geometries (e.g., river lines) to snap to
    tolerance : float, default 1000.0
        Snapping tolerance (in meters, using CRS 3857). 
        Points beyond this distance won't be snapped.
    
    Returns
    -------
    gpd.GeoDataFrame
        Stations with updated geometry and additional columns:
        - 'snapped': boolean flag indicating if point was moved
        - 'snap_distance': distance the point was moved (in meters)
    """
    from shapely.ops import snap
    
    target_crs = 3857
    
    points = gdf_stations.to_crs(target_crs).copy()
    targets = gdf_to_snap.to_crs(target_crs)
    
    targets_union = targets.unary_union
    
    # Store original geometries to calculate distances
    original_geoms = points.geometry.copy()
    
    # Snap each point to the target geometries
    snapped_geoms = [snap(point, targets_union, tolerance) for point in points.geometry]
    
    # Calculate snap distances
    snap_distances = [
        orig.distance(snapped) 
        for orig, snapped in zip(original_geoms, snapped_geoms)
    ]
    
    points['snap_distance'] = snap_distances
    points['geometry'] = snapped_geoms
    points['snapped'] = points['snap_distance'] > 0.01  # Threshold of 1cm for floating point precision
    
    num_snapped = points['snapped'].sum()
    
    print(f'Stations snapped: {num_snapped} out of {len(points)}')

    return points


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
    output_dir: Output directory (defaults to ./outputs)
    file_name: Output filename without extension
    """

    # Setup output directory
    if output_dir is None:
        output_dir = os.path.join(os.getcwd(), 'outputs')
        os.makedirs(output_dir, exist_ok=True)
        
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