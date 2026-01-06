def test_dem_reprojection(dem_array):
    dem_3857 = dem_array.rio.reproject("EPSG:3857")

    assert dem_3857.rio.crs.to_string() == "EPSG:3857"
    assert dem_3857.shape == dem_array.shape
