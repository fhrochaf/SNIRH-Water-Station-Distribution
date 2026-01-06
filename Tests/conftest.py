import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(PROJECT_ROOT))
import pytest
import geopandas as gpd
from shapely.geometry import Point, Polygon, LineString
import numpy as np
import xarray as xr
import rioxarray
import sys
from pathlib import Path




@pytest.fixture
def subbasin_gdf():
    return gpd.GeoDataFrame(
        {"DNS_NU_SUB": [1]},
        geometry=[Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])],
        crs="EPSG:4326"
    )


@pytest.fixture
def stations_gdf():
    return gpd.GeoDataFrame(
        {
            "codigoestacao": [100, 100, 200],
            "Municipio_Nome": ["A", "A", "B"],
            "geometry": [Point(1, 1), Point(1, 1), Point(5, 5)]
        },
        crs="EPSG:4326"
    )


@pytest.fixture
def dem_array():
    data = np.random.rand(20, 20)

    x = np.linspace(0, 10, 20)
    y = np.linspace(0, 10, 20)

    da = xr.DataArray(
        data,
        dims=("y", "x"),
        coords={"x": x, "y": y},
        name="elevation"
    )

    # Register spatial dimensions (CRITICAL)
    da = da.rio.set_spatial_dims(x_dim="x", y_dim="y")

    # Assign CRS
    da = da.rio.write_crs("EPSG:4326")

    return da



@pytest.fixture
def streams_gdf():
    return gpd.GeoDataFrame(
        {
            "strord": [6, 7],
            "geometry": [
                LineString([(1, 1), (4, 4)]),
                LineString([(4, 4), (8, 8)])
            ]
        },
        crs="EPSG:3857"
    )


@pytest.fixture
def catchments_gdf():
    return gpd.GeoDataFrame(
        {
            "value": [1],
            "geometry": [Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])]
        },
        crs="EPSG:3857"
    )
