# some datasets require some geoprocessing (e.g. split, buffer, etc) prior
# to inclusion in the Marxan analysis.

# Each dataset requires different processing. I could probably write some
# common functions, but it is hard to anticipate each case, so this might be
# a bit more verbose than necessary.

import arcpy
import pandas as pd
import numpy as np

root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\spatial'
dir_in = os.path.join(root, '01_original')
gdb_out = os.path.join(root, r'02_processed\dst_processed.gdb')
grid = os.path.join(root, r'03_working\dst_grid.gdb\dst_grid1km')
coastline = os.path.join(root, r'03_working\dst_grid.gdb\dst_coastline')




#######################################
#######################################
# Critical habitat
# extract data for marine species

gdbs = os.path.join(dir_in, 'BC_CriticalHabitat_CB_HabitatEssentiel.gdb')

# check if it intersects with the coastline
arcpy.env.workspace = gdbs
for fc in arcpy.ListFeatureClasses():
    arcpy.Intersect_analysis([fc, coastline], "memory/output")
    result = arcpy.GetCount_management("memory/output")
    count = int(result.getOutput(0))
    if count > 0:
        arcpy.MultipartToSinglepart_management('memory/output', 'memory/multisingle')
        species = fc.split('_', 2)[-1]
        out_name = os.path.join(gdb_out, f'eco_areas_crithab_{species}')
        arcpy.CopyFeatures_management('memory/multisingle', out_name)
    arcpy.Delete_management("memory")

# Delete any features that are just slivers.
arcpy.env.workspace = gdb_out
for fc in arcpy.ListFeatureClasses(wild_card='eco_areas_crithab*'):
    with arcpy.da.UpdateCursor(fc, ['Shape_Area']) as cursor:
        for row in cursor:
            if row[0] < 300:
                cursor.deleteRow()

# Delete any empty fcs
for fc in arcpy.ListFeatureClasses(wild_card='eco_areas_crithab*'):
    result = arcpy.GetCount_management(fc)
    count = int(result.getOutput(0))
    if count == 0:
        arcpy.Delete_management(fc)



#######################################
#######################################
# Substrate model
# Create separate feature classes for mixed, soft, hard
# There is a 20m resolution nearshore dataset and a 100m resolutin deep dataset.
# The deep dataset can be used to fill in gaps of the nearshore dataset.


substrate_near_package = os.path.join(root, r'01_original\substrate_20m_datapackage')
substrate_deep = os.path.join(root, r'01_original\rfsubstrate_100m\Ranger_RF_Substrate_100m.tif')

# convert nearshore to polygons
arcpy.env.workspace = substrate_near_package
for rast in arcpy.ListRasters():
    arcpy.RasterToPolygon_conversion(
        rast,
        os.path.join(gdb_out, 'tmp' + rast.split('.')[0]),
        'NO_SIMPLIFY',
        'SUBSTRATE',
        'SINGLE_OUTER_PART'
    )

# convert deep to polygons
arcpy.env.workspace = gdb_out
arcpy.RasterToPolygon_conversion(
    substrate_deep,
    'tmp_substrate_100m',
    'NO_SIMPLIFY',
    'SUBSTRATE',
    'SINGLE_OUTER_PART'
)

# the boundaries of each region overlap and are not the same
# erase HG with NC
# erase NC with QCS
# erase QCS with SoG
# erase WCVI with SoG
# erase WCVI wth QCS
arcpy.Erase_analysis('tmp_HG_substrate_20m', 'tmp_NCC_substrate_20m', 'tmp_HG_erase')
arcpy.Erase_analysis('tmp_NCC_substrate_20m', 'tmp_QCS_substrate_20m', 'tmp_NCC_erase')
arcpy.Erase_analysis('tmp_QCS_substrate_20m', 'tmp_SOG_substrate_20m', 'tmp_QCS_erase')
arcpy.Erase_analysis('tmp_WCVI_substrate_20m', 'tmp_QCS_substrate_20m', 'tmp_WCVI_erase')
arcpy.Erase_analysis('tmp_SOG_substrate_20m', 'tmp_WCVI_substrate_20m', 'tmp_SOG_erase')

# combine nearshore
fcs = arcpy.ListFeatureClasses(wild_card='*_erase')
arcpy.Merge_management(fcs, 'tmp_substrate_20m')

# erase deep with nearshore
arcpy.Erase_analysis('tmp_substrate_100m', 'tmp_substrate_20m', 'tmp_substrate_100m_erase')

# merge
arcpy.Merge_management(['tmp_substrate_20m', 'tmp_substrate_100m_erase'], 'tmp_substrate_merge')

# clip to grid
arcpy.Clip_analysis('tmp_substrate_merge', grid, 'tmp_substrate_clip')

# rename with soft/mixed/hard
arcpy.CopyFeatures_management('tmp_substrate_clip', 'tmp_substrate_copy')
arcpy.AddField_management('tmp_substrate_copy', 'SUBSTRATE_update', 'TEXT')
with arcpy.da.UpdateCursor('tmp_substrate_copy', ['SUBSTRATE', 'SUBSTRATE_update']) as cursor:
    for row in cursor:
        if row[0] == 'Rock':
            row[1] = 'hard'
        elif row[0] == 'Mixed':
            row[1] = 'mixed'
        elif row[0] == 'Sand' or row[0] == 'Mud':
            row[1] = 'soft'
        cursor.updateRow(row) 

# dissolve by substrate
arcpy.Dissolve_management('tmp_substrate_copy', 'tmp_substrate_dissolve', 'SUBSTRATE_update')
# multipart to singlepart
arcpy.MultipartToSinglepart_management('tmp_substrate_dissolve', 'tmp_substrate_multi')

# split by attribute
arcpy.SplitByAttributes_analysis('tmp_substrate_multi', gdb_out, 'SUBSTRATE_update')

# rename
for fc in ['hard', 'mixed', 'soft']:
    outname = f'eco_phys_substrate_{fc}'
    arcpy.CopyFeatures_management(fc, outname)

# delete intermediate polys
to_delete = arcpy.ListFeatureClasses('tmp*')
for fc in to_delete:
    arcpy.Delete_management(fc)
for fc in ['hard', 'mixed', 'soft']:
    arcpy.Delete_management(fc)



#######################################
#######################################
# Herring spawn index


csv = os.path.join(root, r'01_original\herring_spawn_index\Pacific_herring_spawn_index_data_EN.csv')
arcpy.env.workspace = gdb_out

# remove any with missing coordinates
df = pd.read_csv(csv)
df = df[~df.Latitude.isnull()]

# convert to table
x = np.array(np.rec.fromrecords(df.values))
names = df.dtypes.index.tolist()
x.dtype.names = tuple(names)
arcpy.da.NumPyArrayToTable(x, os.path.join(arcpy.env.workspace, f'temp_table')) # for some reason, you always have to give it the full path

# convert to points
arcpy.XYTableToPoint_management(
    'temp_table',
    'temp_points',
    'Longitude',
    'Latitude',
    coordinate_system = 4326
)

# buffer to create polys
arcpy.env.outputCoordinateSystem = 3005
arcpy.Buffer_analysis(
    'temp_points',
    'temp_buffer',
    100
)

# clip to grid
# (not to coastline. The locations look fairly general, and we may lose some 
# observations by clipping to the coastline)
arcpy.Clip_analysis(
    'temp_buffer',
    grid,
    'eco_fish_herringspawnindex'
)

# Delete
for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)
for fc in arcpy.ListTables('temp*'):
    arcpy.Delete_management(fc)
