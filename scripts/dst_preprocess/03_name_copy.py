
# process/copy datasets to gdb for marxan/prioritizr analysis

import arcpy
import pandas as pd

# download naming scheme from google sheet before
csv = r"C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\scripts\dst_preprocess\DSTpilot_spatialData - naming_scheme.csv"
root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\spatial'
grid = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\spatial\03_working\dst_grid.gdb\dst_grid1km'

names = pd.read_csv(csv)

arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)

for index, row in names.iterrows():

    fc_in = os.path.join(root, row.orig_root, row.orig_folder, row.orig_file)
    fc_out = os.path.join(root, row.processed_root, row.processed_file)

    if arcpy.Exists(fc_out):
        continue

    print(f'Processing {row.processed_file}')

    # clip to grid
    arcpy.Clip_analysis(fc_in, grid, 'memory/clip')

    # dissolve
    arcpy.Dissolve_management('memory/clip', fc_out, multi_part='SINGLE_PART')

    arcpy.Delete_management('memory')


# delete if no features intersected with the grid
arcpy.env.workspace = os.path.join(root, '03_working/dst_eco.gdb')
for fc in arcpy.ListFeatureClasses():
    count = arcpy.GetCount_management(fc)
    if int(count[0]) == 0:
        print(f'Deleting {fc}, no features overlapping')
        arcpy.Delete_management(fc)


# check that all are in the correct coordinate system:
arcpy.env.workspace = os.path.join(root, '03_working/dst_eco.gdb')
for fc in arcpy.ListFeatureClasses():
    spatref = arcpy.Describe(fc).spatialReference.name
    if spatref != 'NAD_1983_BC_Environment_Albers':
        print(fc)

