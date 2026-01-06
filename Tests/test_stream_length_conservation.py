def test_stream_length_conservation(streams_gdf):
    total = streams_gdf.geometry.length.sum()
    assert total > 0
