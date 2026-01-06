import warnings
import geopandas as gpd

def test_station_assignment_integrity(stations_gdf, catchments_gdf):
    with warnings.catch_warnings(record=True) as w:
        gpd.sjoin(stations_gdf, catchments_gdf, predicate="within")
        assert any("CRS mismatch" in str(warn.message) for warn in w)
