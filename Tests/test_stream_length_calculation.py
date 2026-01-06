def test_stream_length_calculation(streams_gdf):
    streams_gdf["length_km"] = streams_gdf.geometry.length / 1000

    assert "length_km" in streams_gdf.columns
    assert streams_gdf["length_km"].gt(0).all()
