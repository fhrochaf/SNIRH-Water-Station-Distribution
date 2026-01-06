def test_dem_produces_catchments(monkeypatch, dem_array, subbasin_gdf):
    """
    Integration safety test:
    Ensures hydrological pipeline can be invoked non-interactively.
    """

    # Mock user input for minimum area
    monkeypatch.setattr("builtins.input", lambda _: "4")

    from catchments.catchments_processing import pyflwdir_subbasins_minarea

    catchments, streams = pyflwdir_subbasins_minarea(
        dem_array.rio.reproject("EPSG:3857"),
        subbasin_gdf,
        min_sto=4,
        save_data=False
    )

    # Weak assertions by design (integration test)
    assert catchments is not None
    assert streams is not None
