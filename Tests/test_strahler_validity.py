def test_strahler_hierarchy(streams_gdf):
    assert streams_gdf["strord"].min() >= 1
    assert streams_gdf["strord"].max() >= streams_gdf["strord"].min()
