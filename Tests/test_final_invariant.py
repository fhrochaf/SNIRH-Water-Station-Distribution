def test_final_analysis_invariant(catchments_gdf):
    gdf = catchments_gdf.copy()
    gdf["count_st_cathcm"] = [3]
    gdf["area_km2"] = [50.0]

    density = gdf["count_st_cathcm"] / gdf["area_km2"]

    assert (density > 0).all()
