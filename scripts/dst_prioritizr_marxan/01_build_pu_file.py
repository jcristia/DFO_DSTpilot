# Calculate area of overlap of each feature with each grid cell
# This will be the input into Prioritizr

import arcpy
import pandas as pd
import numpy as np

root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\spatial'
gdb_features = os.path.join(root, '03_working/dst_eco.gdb')
gdb_grid = os.path.join(root, '03_working/dst_grid.gdb')
new_grid = 'dst_grid_TEMPTEST_PUVSP'

# copy dst_grid_TEMPTEST
# add uID field and calculate as OBJECTID
# recalculate AREA and COST fields to be 1000000
arcpy.env.workspace = gdb_grid
arcpy.CopyFeatures_management('dst_grid_TEMPTEST', new_grid)
arcpy.AddField_management(new_grid, 'uID', 'LONG')
arcpy.CalculateField_management(new_grid, 'uID', '!OBJECTID!')
arcpy.CalculateField_management(new_grid, 'AREA', 1000000)
arcpy.CalculateField_management(new_grid, 'COST', 1000000)
arcpy.DeleteField_management(new_grid, 'UNIT_ID')


# use Tabulate Intersection tool to assign area of overlapping feature to each 
# grid cell.
cursor = arcpy.da.SearchCursor(new_grid, ['uID'])
df_all = pd.DataFrame(data=[row for row in cursor], columns=['uID'])
arcpy.env.workspace = gdb_features
for fc in arcpy.ListFeatureClasses():
    print(f'Processing {fc}')

    arcpy.TabulateIntersection_analysis(
        os.path.join(gdb_grid, new_grid),
        ['uID'],
        fc,
        'temp_table',
    )

    field_names = ['uID', 'AREA']
    cursor = arcpy.da.SearchCursor('temp_table', field_names)
    df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    df = df.rename(columns={'AREA':fc})
    df_all = df_all.merge(df, how='left', on='uID')
    arcpy.Delete_management('temp_table')

# join df back to grid
x = np.array(np.rec.fromrecords(df_all.values))
names = df_all.dtypes.index.tolist()
x.dtype.names = tuple(names)
arcpy.da.NumPyArrayToTable(x, os.path.join(gdb_grid, f'temp_tbl'))
arcpy.env.qualifiedFieldNames = False
arcpy.env.workspace = gdb_grid
jg = arcpy.AddJoin_management(new_grid, 'uID', 'temp_tbl', 'uID')
arcpy.CopyFeatures_management(jg, f'{new_grid}_joined')
arcpy.Delete_management('temp_tbl')
arcpy.DeleteField_management(f'{new_grid}_joined', ['OBJECTID_1', 'UID_1'])
arcpy.Delete_management(new_grid)


# Other option: intersect, dissolve, to pandas table, but would require more cleanup.
# I also tried the Apportion tool but ran into field length limit errors.
