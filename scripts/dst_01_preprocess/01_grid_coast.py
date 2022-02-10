# create 1km grid from existing 2km grid


import arcpy


base = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\spatial'
coastline = os.path.join(base, r'01_original\Pacific_fine_scale_ocean_fromTRIM_20210615.gdb\Pacific_Ocean_fromTRIM_20210609')
grid_2km = os.path.join(base, r'01_original\DSTpilot_grids.gdb\GRID_MARX_planningUnits2x2_bc')
gdb = os.path.join(base, r'03_working\dst_grid.gdb')
ss_bioregion_EXTENT = 'DFO_Marine_Bioregions_SouthernShelf_CUSTOMEXTENT'
arcpy.env.workspace = gdb


# copy coastline for easy reference in working folder
arcpy.CopyFeatures_management(coastline, 'dst_coastline')


# split 2km grid into 1km grid
arcpy.SubdividePolygon_management(
    grid_2km, 
    'dst_grid1km_01divide', 
    method = 'NUMBER_OF_EQUAL_PARTS', 
    num_areas = 4,
    split_angle = 0,
    subdivision_type='STACKED_BLOCKS'
    )

# clip to Southern Shelf Bioregion
# HOWEVER, the coastline used for bioregions is not as high resolution as the
# TRIM coastline. If I clip with just this then some PUs will be removed that
# should not be. Therefore, I need to get a general extent of the SSBR to clip
# to.
# I MANUALLY drew this, by just making a general shape that includes all of the
# coastline and then traced the pelagic outline of the SSBR.
arcpy.Clip_analysis('dst_grid1km_01divide', ss_bioregion_EXTENT, 'dst_grid1km_02clip')


# remove any grid cells that don't overlap with the coastline
# (this takes a few minutes to run, but it does work)
arcpy.MakeFeatureLayer_management('dst_grid1km_02clip', 'temp_lyr')
sel_grid = arcpy.SelectLayerByLocation_management(
    'temp_lyr',
    'INTERSECT',
    coastline,
    selection_type='NEW_SELECTION'
)
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)
arcpy.CopyFeatures_management(sel_grid, 'dst_grid1km')

arcpy.Delete_management('temp_lyr')
arcpy.Delete_management('dst_grid1km_01divide')
arcpy.Delete_management('dst_grid1km_02clip')