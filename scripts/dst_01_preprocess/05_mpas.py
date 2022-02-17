# quantify overlap of MPAs/RCAs and planning units


import arcpy

root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot\spatial'
dir_in = os.path.join(root, '01_original')
gdb_out = os.path.join(root, r'03_working/dst_mpas.gdb')
grid = os.path.join(root, r'03_working\dst_grid.gdb\dst_grid1km_PUVSP')
mpas = os.path.join(dir_in, r'CPCAD-BDCAPC_Dec2021.gdb\CPCAD_BDCAPC_Dec2021')
rcas = os.path.join(dir_in, r'RCA\rcas2005.shp')
arcpy.env.workspace = gdb_out



### MPAs
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)
arcpy.PairwiseClip_analysis(mpas, grid, 'temp_mpas_clip')
# get unique list of MPA names
ms = []
with arcpy.da.SearchCursor('temp_mpas_clip', ['NAME_E']) as cursor:
    for row in cursor:
        if row[0] not in ms:
            ms.append(row[0])
# We will only include terrestrial areas if they are adjacent to a marine component.
# Identify areas that are not:
terrestrial_remove = []
for m in ms:
    where = f"""NAME_E = '{m}' And BIOME = 'M'"""
    arcpy.MakeFeatureLayer_management('temp_mpas_clip', 'temp_out', where_clause=where)
    count = arcpy.GetCount_management('temp_out')
    if int(count[0]) == 0:
        print(f'{m} no selection')
        terrestrial_remove.append(m)
    arcpy.Delete_management('temp_out')
# So there are actually a few marine parks that look like they were mistakenly
# coded as 'T'. So a simple select will not work for these.
# Also, once I scan some of these on the map, I think we can be ok with just
# clipping to the grid and keeping any terrestrial parks selected. It seems like
# most terrestrial areas selected are adjacent to a marine reserve, so it is ok
# for that PU to be locked in.
# However, I think Carrie should review this.
# For now, just copy a final version.
arcpy.CopyFeatures_management('temp_mpas_clip', 'mpas_polygons')
arcpy.Delete_management('temp_mpas_clip')



### RCAs
arcpy.Clip_analysis(rcas, grid, 'rcas_polygons')



### Merge MPAs and RCAs
arcpy.Merge_management(['mpas_polygons', 'rcas_polygons'], 'mpas_rcas_polygons')
for field in arcpy.ListFields('mpas_rcas_polygons'):
    if not field.required and field.name not in ['NAME_E', 'NAME']:
        arcpy.DeleteField_management('mpas_rcas_polygons', field.name)
arcpy.AddField_management('mpas_rcas_polygons', 'name_mpa_rca', 'TEXT', field_length=200)
arcpy.AddField_management('mpas_rcas_polygons', 'mpa_rca', 'TEXT')
with arcpy.da.UpdateCursor('mpas_rcas_polygons', ['NAME_E', 'NAME', 'name_mpa_rca', 'mpa_rca']) as cursor:
    for row in cursor:
        if row[0] == None:
            row[2] = row[1]
            row[3] = 'rca'
        else:
            row[2] = row[0]
            row[3] = 'mpa'
        cursor.updateRow(row)
arcpy.DeleteField_management('mpas_rcas_polygons', 'NAME_E')
arcpy.DeleteField_management('mpas_rcas_polygons', 'NAME')



### PU overlap for merged dataset
arcpy.Intersect_analysis(['mpas_rcas_polygons', grid], 'temp_mpas_rcas_marxan_intersect')
arcpy.Dissolve_management('temp_mpas_rcas_marxan_intersect', 'mpas_rcas_marxan', ['uID'])
arcpy.Delete_management('temp_mpas_rcas_marxan_intersect')



### Threshold
# when deciding to lock in a planning unit, perhpas we set a threshold of
# overlap so that a sliver doesn't cause a whole unit to be locked in.
# We could do this in R, but do a version of it here since we will need it for
# Marxan.
arcpy.MakeFeatureLayer_management(
    'mpas_rcas_marxan',
    'temp_lyr',
    'Shape_Area > 150000'
)
arcpy.CopyFeatures_management('temp_lyr', 'mpas_rcas_marxan_threshold150k')