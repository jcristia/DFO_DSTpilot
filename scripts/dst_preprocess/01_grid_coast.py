# create 1km grid from existing 2km grid


import arcpy


base = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\spatial'
coastline = os.path.join(base, r'01_original\Pacific_fine_scale_ocean_fromTRIM_20210615.gdb\Pacific_Ocean_fromTRIM_20210609')
grid_2km = os.path.join(base, r'01_original\DSTpilot_grids.gdb\GRID_MARX_planningUnits2x2_bc')
gdb = os.path.join(base, r'03_working\dst_grid.gdb')
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


# remove any grid cells that don't overlap with the coastline
# (this takes a few minutes to run, but it does work)
arcpy.MakeFeatureLayer_management('dst_grid1km_01divide', 'temp_lyr')
sel_grid = arcpy.SelectLayerByLocation_management(
    'temp_lyr',
    'INTERSECT',
    coastline,
    selection_type='NEW_SELECTION'
)
arcpy.CopyFeatures_management(sel_grid, 'dst_grid1km')
arcpy.Delete_management('temp_lyr')
arcpy.Delete_management('dst_grid1km_01divide')