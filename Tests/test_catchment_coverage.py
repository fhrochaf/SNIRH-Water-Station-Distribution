def test_catchments_cover_subbasin(catchments_gdf, subbasin_gdf):
    union_catchments = catchments_gdf.unary_union
    subbasin_geom = subbasin_gdf.geometry.iloc[0]

    assert union_catchments.contains(subbasin_geom.centroid)
