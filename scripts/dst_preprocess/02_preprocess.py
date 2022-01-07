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
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)
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



#######################################
#######################################
# Harbour seal haulouts
# buffer following BCMCA methodology

csv = os.path.join(root, r'01_original\harbour_seal\Harbour_Seal_Counts_from_Strait_of_Georgia_Haulout_Locations.csv')
# some points have a positive longitude
pts = pd.read_csv(csv)
pts['Longitude'] = pts.Longitude.abs() * -1
pts.to_csv(os.path.join(root, r'01_original\harbour_seal\edited.csv'))
csv = os.path.join(root, r'01_original\harbour_seal\edited.csv')

arcpy.env.workspace = gdb_out
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)
arcpy.XYTableToPoint_management(
    csv,
    'temp_points',
    'Longitude',
    'Latitude',
    coordinate_system = 4326
)

os.remove(csv)

arcpy.Buffer_analysis('temp_points', 'temp_buffer', 200)
arcpy.Dissolve_management(
    'temp_buffer', 
    'eco_mammals_harboursealhaulouts', 
    multi_part='SINGLE_PART')

arcpy.Delete_management('temp_points')
arcpy.Delete_management('temp_buffer')



#######################################
#######################################
# Steller sea lion haulouts
# Rookeries and haulouts into separate datasets
# Buffer following BCMCA methodology
# 200m haulouts, 15km rookeries

csv = os.path.join(root, r'01_original\steller_sealion\Steller_Sea_Lion_Summer_counts_from_Haulout_Locations.csv')

# from metadata:
# Y: Year-round haulout site
# W: Winter haulout site
# A: Breeding rookeries

# other codes in dataset: Null, ?, R, W/Y, Y/A, Y/R, Y/W
# also, there are some records without coordinates

df = pd.read_csv(csv)
haul_rook = [['Y', 'W', 'W/Y', 'Y/W'], ['R', 'Y/A', 'Y/R']]
names = ['haulouts', 'rookeries']
buffs = [200, 15000]
arcpy.env.workspace = gdb_out
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)

for hr, name, buff in zip(haul_rook, names, buffs):
    
    df_temp = df[df['SITE TYPE'].isin(hr) & ~df.LATITUDE.isnull()]
    x = np.array(np.rec.fromrecords(df_temp.values))
    names = df_temp.dtypes.index.tolist()
    x.dtype.names = tuple(names)
    arcpy.da.NumPyArrayToTable(x, os.path.join(arcpy.env.workspace, f'temp_table'))

    arcpy.XYTableToPoint_management(
        'temp_table',
        'temp_points',
        'LONGITUDE',
        'LATITUDE',
        coordinate_system = 4326
    )

    arcpy.Buffer_analysis('temp_points', 'temp_buffer', buff)
    arcpy.Dissolve_management(
        'temp_buffer', 
        f'eco_mammals_stellersealion{name}', 
        multi_part='SINGLE_PART')

    arcpy.Delete_management('temp_points')
    arcpy.Delete_management('temp_table')
    arcpy.Delete_management('temp_buffer')



#######################################
#######################################
# Fjords
# Select from oceanograhic regions dataset

arcpy.env.workspace = gdb_out
fc = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\mpatt_eco_coarse_oceanographicRegions_data')
where = "Id = 11"
arcpy.MakeFeatureLayer_management(
    fc, 
    'temp_lyr', 
    where
    )
arcpy.CopyFeatures_management('temp_lyr', 'eco_phys_fjords')
arcpy.Delete_management('temp_lyr')



#######################################
#######################################
# Kelp biobands
# buffer, dissolve, merge, clip

canopy = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\bc_ShoreZone_canopy_kelps_march2020')
under = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\bc_ShoreZone_under_kelp_march2020')
arcpy.env.workspace = gdb_out

# Arc might crash on these. Just keep rerunning.

arcpy.Buffer_analysis(canopy, 'temp_buff_can', 50, line_end_type='FLAT')
arcpy.Buffer_analysis(under, 'temp_buff_und', 50, line_end_type='FLAT')
arcpy.Merge_management(['temp_buff_can', 'temp_buff_und'], 'temp_merge')
arcpy.Dissolve_management('temp_merge', 'temp_dissolve', multi_part='SINGLE_PART')
arcpy.Clip_analysis('temp_dissolve', coastline, 'eco_plants_kelpbiobands')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)



#######################################
#######################################
# Eelgrass polygons and biobands
# buffer biobands, merge, dissolve and singlepart, clip to coastline

eg_poly = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\eelgrass_BC_polygons_explode')
eg_line = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\bc_ShoreZone_eelgrass_march2020')
arcpy.env.workspace = gdb_out

# polygons have significant errors in them
arcpy.CopyFeatures_management(eg_poly, 'temp_copy')
arcpy.RepairGeometry_management('temp_copy')
arcpy.CheckGeometry_management('temp_copy')

arcpy.Buffer_analysis(eg_line, 'temp_buff', 50, line_end_type='FLAT')
arcpy.Merge_management(['temp_buff', 'temp_copy'], 'temp_merge')
arcpy.PairwiseDissolve_analysis('temp_merge', 'temp_dissolve', multi_part='SINGLE_PART') # crashed if using normal dissolve
arcpy.PairwiseClip_analysis('temp_dissolve', coastline, 'eco_plants_eelgrass')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)
arcpy.Delete_management('temp_copy_CheckGeometry')


#######################################
#######################################
# Surfgrass biobands
# buffer

sg_line = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\mpatt_eco_plants_surfgrass_biobands_data')
arcpy.env.workspace = gdb_out

arcpy.Buffer_analysis(sg_line, 'temp_buff', 25, line_end_type='FLAT')
arcpy.Dissolve_management('temp_buff', 'temp_dissolve', multi_part='SINGLE_PART')
arcpy.Clip_analysis('temp_dissolve', coastline, 'eco_plants_surfgrass')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)




#######################################
#######################################
# High rugosity
# Extract out high values
# BP: BCMCA calculated the natural logarithm of the rugosity layer, and then 
# classified rugosity into quantiles, and the top 20% of the values were 
# identified as "high rugosity".
# https://bcmca.ca/datafiles/individualfiles/bcmca_eco_physical_highrugosity_atlas.pdf


# STILL NEED TO CLIP TO GRID AFTER MERGE, BEFORE I TAKE TOP 20%. THIS MIGHT BE
# ABLE TO TAKE THE PLACE OF SPECIFYING A SUBSET WHEN GOING TO NUMPY.
# IF I CLIP TO A RECTANGLE, THEN I WILL ALSO NEED TO ASSIGN NODATA VALUES TO
# AREAS LIKE WCVI AND RIVERS THAT SHOULD NOT BE INCLUDED.
# If I end up needing to calculate the percentile over the entire grid, then
# I will need to break it up into individual raster and bring them into numpy:
# https://community.esri.com/t5/python-questions/how-to-calculate-the-percentile-for-each-cell-from/td-p/130653
# I just need the 80% value, so after that, I can then just apply it to the one
# study area raster.


rugs = os.path.join(dir_in, 'DST_PilotEcologicalData_BP.gdb')
arcpy.env.workspace = rugs

# combine regional rasters
from arcpy.ia import *
rast1 = arcpy.Raster('rugosity_shelfsalishsea')
rast2 = arcpy.Raster('rugosity_nearshore_WCVI')
rast3 = arcpy.Raster('rugosity_nearshore_qcs')
rast4 = arcpy.Raster('rugosity_nearshore_ncc')
rast5 = arcpy.Raster('rugosity_nearshore_hg')
rast_reg = arcpy.ia.Merge([rast1, rast2, rast3, rast4, rast5], 'MEAN')

# resample nsb_ssb raster to 20m cell size and snap to others
arcpy.env.snapRaster = 'rugosity_shelfsalishsea'
arcpy.Resample_management(
    'rugosity_nsb_ssb',
    os.path.join(gdb_out, 'temp_resample'),
    cell_size = 20,
    resampling_type='NEAREST'
)

# combine, where they overlap take the regional nearshore raster value
arcpy.env.workspace = gdb_out
rast6 = arcpy.Raster('temp_resample')
rast_all = arcpy.ia.Merge([rast_reg, rast6], 'FIRST')
rast_all.save('temp_merge')

# It looks like the only way to properly calculate percentiles is to use numpy. 
# However, there is a size limit, which this raster exceeds.
# Therefore, specify a subset.

rast_all = arcpy.Raster('temp_merge')
checkymin = rast_all.extent.YMin
checkxmax = rast_all.extent.XMax

# bottom left corner
bottom = 347000
left = 1033000
lowerleft = arcpy.Point(left, bottom)

# top right corner
top = 570000
right = 1240000

# RasterToNumPyArray requires specifiy number of columns and rows instead of
# an upper right corner. This does't work out perfectly though. Double check.
ncols = round((top-bottom)/20)
nrows = round((right-left)/20)

arr = arcpy.RasterToNumPyArray(rast_all, lowerleft, ncols, nrows)
# convert 0 values to nan so that they aren't considered (I was trying to do 
# this with masking before, but the nanpercentile doesn't work properly with
# that, even though it appears that it does).
arr_nan = np.where(arr == 0, np.nan, arr)
n_80 = np.nanpercentile(arr_nan, 80) # this is the 80th percentile value
# now, get just the values greater than that
arr_80 = np.where(arr_nan < n_80, np.nan, arr_nan)
# convert nan to 0
arr_out = np.nan_to_num(arr_80)

# save a version with values
rast_arr = arcpy.NumPyArrayToRaster(
    arr_out,
    lowerleft, 
    20, 20, 
    value_to_nodata=0)
rast_arr.save('temp_perc_values')

# save a version with just 1/0 to convert to polys
arr_80 = np.where(arr_nan < n_80, np.nan, arr_nan)
arr_out = np.where(arr_80 > 0, 1, 0)
rast_arr = arcpy.NumPyArrayToRaster(
    arr_out,
    lowerleft, 
    20, 20, 
    value_to_nodata=0)
rast_arr.save('temp_perc_novalues')

# convert to polys
arcpy.RasterToPolygon_conversion(
    'temp_perc_novalues',
    'eco_phys_highrugosity',
    simplify='NO_SIMPLIFY',
)

for rast in arcpy.ListRasters('temp*'):
    arcpy.Delete_management(rast)




#######################################
#######################################
# Cold seeps
# buffer points by 500 meters


cs_points = os.path.join(dir_in, 'DSTpilot_ecologicalData.gdb/Confirmed_cold_seeps_DFO_2018')
arcpy.env.workspace = gdb_out

arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)
arcpy.Buffer_analysis(
    cs_points,
    'memory/temp_buff',
    500,
    dissolve_option='NONE'
)
arcpy.Dissolve_management(
    'memory/temp_buff',
    'eco_phys_coldseeps',
    multi_part='SINGLE_PART'
)
arcpy.Delete_management('memory')




#######################################
#######################################
# Geomorphic units
# split by attribute

geomorph_units = os.path.join(dir_in, 'DST_PilotEcologicalData_BP.gdb/BC_Geomorphic_Units')
arcpy.env.workspace = gdb_out

# view list of unique values
geomorph = []
geogreene = []
with arcpy.da.SearchCursor(geomorph_units, ['Geomorph', 'GeoGreene']) as cursor:
    for row in cursor:
        if row[0] not in geomorph:
            geomorph.append(row[0])
        if row[1] not in geogreene:
            geogreene.append(row[1])

# do just for slop and shelf for now
atts = {
    'slope_canyonfloor':'Slope, Canyon floor',
    'slope_ridge':'Slope, Ridge',
    'shelf_wallsloping':'Shelf, Wall, sloping',
    'shelf_crest':'Shelf, Crest',
    'shelf_depressionfloor':'Shelf, Depression floor',
    'shelf_mound':'Shelf, Mound',
    'shelf_depression':'Shelf, Depression'
}

for key in atts:

    where = f""""GeoGreene" = '{atts[key]}'"""
    arcpy.MakeFeatureLayer_management(
        geomorph_units,
        'temp_out',
        where_clause=where
    )
    out_name = f'eco_phys_geomorph_{key}'
    arcpy.CopyFeatures_management('temp_out', out_name)
    arcpy.Delete_management('temp_out')




#######################################
#######################################
# Marine bird colonies

gdb_birds = os.path.join(dir_in, 'mpatt_eco_birds_colonies.gdb')
arcpy.env.workspace = gdb_birds
fcs = arcpy.ListFeatureClasses()

for fc in fcs:
    bird = fc.split('_')[3]
    out_name = os.path.join(gdb_out, f'eco_birds_{bird}_colonies')
    arcpy.Dissolve_management(fc, out_name, multi_part='SINGLE_PART')




#######################################
#######################################
# Important areas
# Each group requires difference geoprocessing so must be done separately

ia_dir = os.path.join(dir_in, 'ImportantAreas_Databases')
arcpy.env.workspace = gdb_out

## Birds ##
walk = arcpy.da.Walk(ia_dir, datatype="FeatureClass", type="Polygon")
fcs = []
for dirpath, dirnames, filenames in walk:
    for filename in filenames:
        if 'birds' in filename:
            fcs.append(os.path.join(dirpath, filename))
arcpy.Merge_management(fcs, 'memory/merge')
arcpy.Dissolve_management('memory/merge', 'eco_areas_ia_birds', multi_part='SINGLE_PART')
arcpy.Delete_management('memory')


## Cetaceans ##
# multiple species across 3 geodatabases
# I manually changed names across the geodatabases so that I can easily match
# them. This saves a lot of code.
walk = arcpy.da.Walk(ia_dir, datatype="FeatureClass", type="Polygon")
fcs = []
for dirpath, dirnames, filenames in walk:
    for filename in filenames:
        if 'cetaceans' in dirpath:
            fcs.append(os.path.join(dirpath, filename))
species = [] # get unique list
for fc in fcs:
    species.append(fc.split('_')[-2])
species_set = set(species)
species = list(species_set)

for spec in species:
    merge_list = []
    for fc in fcs:
        if spec in fc:
            merge_list.append(fc)
    arcpy.Merge_management(merge_list, 'memory/merge')
    arcpy.Dissolve_management('memory/merge', f'eco_cetaceans_ia_{spec}', multi_part='SINGLE_PART')
    arcpy.Delete_management('memory')


## coralsponge ##
# species are all mixed
# overlaps with other sponge reef dataset. May want to erase with this.
# just merge for now until I know what the team wants to do
walk = arcpy.da.Walk(ia_dir, datatype="FeatureClass", type="Polygon")
fcs = []
for dirpath, dirnames, filenames in walk:
    for filename in filenames:
        if 'coralsponge' in dirpath:
            fcs.append(os.path.join(dirpath, filename))
arcpy.Merge_management(fcs, 'memory/merge')
arcpy.Dissolve_management('memory/merge', f'eco_inverts_ia_coralsponge', multi_part='SINGLE_PART')
arcpy.Delete_management('memory')


## fish ##
# changed one name manually (sixgill-ed shark)
walk = arcpy.da.Walk(ia_dir, datatype="FeatureClass", type="Polygon")
fcs = []
for dirpath, dirnames, filenames in walk:
    for filename in filenames:
        if 'fish' in dirpath:
            fcs.append(os.path.join(dirpath, filename))
species = [] # get unique list
for fc in fcs:
    species.append(fc.split('_')[-2])
species_set = set(species)
species = list(species_set)
species = sorted(species)

for spec in species:
    merge_list = []
    for fc in fcs:
        if spec in fc:
            merge_list.append(fc)
    arcpy.Merge_management(merge_list, 'memory/merge')
    arcpy.Dissolve_management('memory/merge', f'eco_fish_ia_{spec}', multi_part='SINGLE_PART')
    arcpy.Delete_management('memory')


## invertebrates ##
# split out bivalves first
arcpy.env.workspace = os.path.join(ia_dir, 'ia_invertebrates_sog.gdb')
bv_sog = 'ia_bivalves_sog'
arcpy.MakeFeatureLayer_management(bv_sog, 'temp', """species = 'Butter Clam'""")
arcpy.CopyFeatures_management('temp', 'ia_butterclam_sog')
arcpy.Delete_management('temp')
arcpy.MakeFeatureLayer_management(bv_sog, 'temp', """species = 'Manilla Clam'""")
arcpy.CopyFeatures_management('temp', 'ia_manilaclam_sog')
arcpy.Delete_management('temp')
arcpy.MakeFeatureLayer_management(bv_sog, 'temp', """species = 'Pacific oyster'""")
arcpy.CopyFeatures_management('temp', 'ia_oyster_sog')
arcpy.Delete_management('temp')

arcpy.env.workspace = os.path.join(ia_dir, 'ia_invertebrates_wcvi.gdb')
bv_wvi = 'ia_bivalves_wcvi'
arcpy.MakeFeatureLayer_management(bv_wvi, 'temp', """species = 'Razor Clam'""")
arcpy.CopyFeatures_management('temp', 'ia_razorclam_wcvi')
arcpy.Delete_management('temp')
arcpy.MakeFeatureLayer_management(bv_wvi, 'temp', """species IN ('Olympia oyster', 'Pacific oyster')""")
arcpy.CopyFeatures_management('temp', 'ia_oyster_wcvi')
arcpy.Delete_management('temp')

# then go through each one, but remove bivalves from list
walk = arcpy.da.Walk(ia_dir, datatype="FeatureClass", type="Polygon")
fcs = []
for dirpath, dirnames, filenames in walk:
    for filename in filenames:
        if 'invertebrates' in dirpath:
            fcs.append(os.path.join(dirpath, filename))
species = [] # get unique list
for fc in fcs:
    species.append(fc.split('_')[-2])
species_set = set(species)
species = list(species_set)
species = sorted(species)
species.remove('bivalves')

arcpy.env.workspace = gdb_out
for spec in species:
    merge_list = []
    for fc in fcs:
        if spec in fc:
            merge_list.append(fc)
    arcpy.Merge_management(merge_list, 'memory/merge')
    arcpy.Dissolve_management('memory/merge', f'eco_inverts_ia_{spec}', multi_part='SINGLE_PART')
    arcpy.Delete_management('memory')

# delete clam/oyster datasets after merge
arcpy.env.workspace = os.path.join(ia_dir, 'ia_invertebrates_sog.gdb')
arcpy.Delete_management('ia_butterclam_sog')
arcpy.Delete_management('ia_manilaclam_sog')
arcpy.Delete_management('ia_oyster_sog')
arcpy.env.workspace = os.path.join(ia_dir, 'ia_invertebrates_wcvi.gdb')
arcpy.Delete_management('ia_razorclam_wcvi')
arcpy.Delete_management('ia_oyster_wcvi')


## other vertebrates ##
# some fc names manually changed to match
walk = arcpy.da.Walk(ia_dir, datatype="FeatureClass", type="Polygon")
fcs = []
for dirpath, dirnames, filenames in walk:
    for filename in filenames:
        if 'other_vert' in dirpath:
            fcs.append(os.path.join(dirpath, filename))
species = [] # get unique list
for fc in fcs:
    species.append(fc.split('_')[-2])
species_set = set(species)
species = list(species_set)
species = sorted(species)

for spec in species:
    merge_list = []
    for fc in fcs:
        if spec in fc:
            merge_list.append(fc)
    arcpy.Merge_management(merge_list, 'memory/merge')
    if spec == 'leatherbackseaturtle':
        outname = f'eco_reptiles_ia_{spec}'
    else:
        outname = f'eco_mammals_ia_{spec}'
    arcpy.Dissolve_management('memory/merge', outname, multi_part='SINGLE_PART')
    arcpy.Delete_management('memory')

