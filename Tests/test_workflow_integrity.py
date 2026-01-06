def test_stream_length_aggregation(streams_gdf, catchments_gdf):
    joined = streams_gdf.sjoin(
        catchments_gdf,
        how="left",
        predicate="within"
    )

    joined["length_km"] = joined.geometry.length / 1000
    total_length = joined.groupby("index_right")["length_km"].sum()

    assert total_length.iloc[0] > 0
