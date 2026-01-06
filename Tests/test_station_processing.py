def test_station_deduplication(stations_gdf):
    dedup = (
        stations_gdf
        .groupby("codigoestacao", as_index=False)
        .last()
    )

    assert len(dedup) == 2
    assert set(dedup["codigoestacao"]) == {100, 200}
