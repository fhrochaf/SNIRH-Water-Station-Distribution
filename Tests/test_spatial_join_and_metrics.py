def test_station_catchment_join(stations_gdf, subbasin_gdf):
    joined = stations_gdf.sjoin(
        subbasin_gdf,
        how="left",
        predicate="within"
    )

    assert joined.notnull().any(axis=None)
