# Calculate area of overlap of each feature with each grid cell
# This will be the input into Prioritizr
# Prioritizr does not like big or small numbers. It causes instability when 
# solving. Therefore, convert everything to a % within each cell. This preserves
# the most values, otherwise we will be losing a lot of smaller values when we
# convert to kilometers squared.

import arcpy
import pandas as pd
import numpy as np

root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\spatial'
gdb_features = os.path.join(root, '03_working/dst_ecohu.gdb')
gdb_grid = os.path.join(root, '03_working/dst_grid.gdb')
new_grid = 'dst_grid1km_PUVSP'



# copy dst_grid_TEMPTEST
# add uID field and calculate as OBJECTID
# recalculate AREA and COST fields to be 100
arcpy.env.workspace = gdb_grid
arcpy.CopyFeatures_management('dst_grid1km', new_grid)
arcpy.AddField_management(new_grid, 'uID', 'LONG')
arcpy.CalculateField_management(new_grid, 'uID', '!OBJECTID!')
arcpy.CalculateField_management(new_grid, 'AREA', 100)
arcpy.CalculateField_management(new_grid, 'COST', 100)
arcpy.DeleteField_management(new_grid, 'UNIT_ID')


# use Tabulate Intersection tool to assign area of overlapping feature to each 
# grid cell. This is just for area based features.
cursor = arcpy.da.SearchCursor(new_grid, ['uID'])
df_all = pd.DataFrame(data=[row for row in cursor], columns=['uID'])
arcpy.env.workspace = gdb_features
fcs = [fc for fc in arcpy.ListFeatureClasses() if not fc.endswith('_d')]
for fc in fcs:

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
    df[fc] = df[fc] / 1000000.0 * 100 # convert to %
    # Deal with small numbers. Round to 4 decimal places.
    # This means that anything less than this will be 0. However, if you view
    # the data, they are just extremely small slivers at this scale.
    df[fc] = df[fc].round(4)
    df_all = df_all.merge(df, how='left', on='uID')
    df_all[fc] = df_all[fc].fillna(0) # replace nan with 0
    arcpy.Delete_management('temp_table')


# deal with non-area based features
for fc in arcpy.ListFeatureClasses('*_d'):
    print(f'Processing {fc}')

    # assign values to grid
    arcpy.env.qualifiedFieldNames = False
    arcpy.SpatialJoin_analysis(
        os.path.join(gdb_grid, new_grid), 
        fc, 
        'temp_table', 
        'JOIN_ONE_TO_MANY', 
        'KEEP_COMMON', 
        match_option='HAVE_THEIR_CENTER_IN'
        )

    for field in arcpy.ListFields('temp_table'):
        if field.name.endswith('_DST'):
            narea_field = field.name
    for field in arcpy.ListFields('temp_table'):
        if not field.required and field.name not in ['uID', narea_field]:
            arcpy.DeleteField_management('temp_table', field.name)

    field_names = ['uID', narea_field]
    cursor = arcpy.da.SearchCursor('temp_table', field_names)
    df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    df = df.rename(columns={narea_field:fc})

    # for non-area features, scale to a max of 100
    # divide by max, multiply by 100, round to 4 decimal places
    df[fc] = df[fc] / df[fc].max() * 100
    df[fc] = df[fc].round(4)
    df_all = df_all.merge(df, how='left', on='uID')
    df_all[fc] = df_all[fc].fillna(0) # replace nan with 0
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

arcpy.Rename_management(f'{new_grid}_joined', new_grid)


# Other option: intersect, dissolve, to pandas table, but would require more cleanup.
# I also tried the Apportion tool but ran into field length limit errors.
