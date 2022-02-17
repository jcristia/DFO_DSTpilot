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
# Critical habitat (ECCC dataset)
# extract data for marine species

gdbs = os.path.join(dir_in, 'BC_CriticalHabitat_CB_HabitatEssentiel.gdb')

# check if it intersects with the grid
arcpy.env.workspace = gdbs
for fc in arcpy.ListFeatureClasses():
    arcpy.Intersect_analysis([fc, grid], "memory/output")
    result = arcpy.GetCount_management("memory/output")
    count = int(result.getOutput(0))
    if count > 0:
        arcpy.MultipartToSinglepart_management('memory/output', 'memory/multisingle')
        species = fc.split('_', 2)[-1]
        out_name = os.path.join(gdb_out, f'crithab_{species}')
        arcpy.CopyFeatures_management('memory/multisingle', out_name)
    arcpy.Delete_management("memory")


# MANUALLY:
# Inspect each feature and look up what kind of species it is.
# Determine if it should be included. If it is terrestrial, it needs to have
# an exclusive coastal association.
to_remove = [
    'Bartramia_stricta', 
    'Carex_tumulicola', 
    'Centaurium_muehlenbergii', 
    'Chrysemys_picta_bellii', 
    'Epilobium_densiflorum', 
    'Lasthenia_glaberrima', 
    'Leptogium_platynum', 
    'Lupinus_densiflorus',
    'Prophysaon_coeruleum', 
    'Ranunculus_californicus',
    'Hypogymnia_heterophylla',
    'Castilleja_victoriae',
    'Limnanthes_macounii',
    'Microseris_bigelovi'
]

# Delete features from list
arcpy.env.workspace = gdb_out
for fc in arcpy.ListFeatureClasses(wild_card='crithab*'):
    for species in to_remove:
        if species in fc:
            arcpy.Delete_management(fc)
            continue

# Rename
spec_type = {
    'birds': ['Brachyramphus_marmoratus'],
    'inverts': ['Anarta_edwardsii'],
    'plants': [
        'Abronia_umbellata',
        'Camissonia_contorta',
        ]
}
common_name = { # not all have common names
    'Brachyramphus_marmoratus':'marbledmurrelet',
    'Anarta_edwardsii': 'edwardsbeachmoth',
    'Abronia_umbellata': 'pinksandverbena',
    'Camissonia_contorta': 'plainseveningprimrose',
}

arcpy.env.workspace = gdb_out
for fc in arcpy.ListFeatureClasses(wild_card='crithab*'):
    species = species = fc.split('_', 1)[-1]
    for t in spec_type:
        if species in spec_type[t]:
            species_type=t
    cname = common_name[species]
    #species = ''.join(species.split('_')).lower()
    outname = f'eco_{species_type}_{cname}'
    arcpy.Rename_management(fc, outname)



#######################################
#######################################
# Critical habitat (DFO dataset)
# Extract out killer whale data. It is the only species overlapping with our
# study area.

gdb = os.path.join(dir_in, 'CriticalHabitat_FGP.gdb')
arcpy.env.workspace = gdb
ch = 'DFO_SARA_CritHab_2021_FGP_EN'
arcpy.MakeFeatureLayer_management(ch, 'temp', """Population_EN = 'Northeast Pacific Southern Resident'""")
arcpy.CopyFeatures_management('temp', os.path.join(gdb_out, 'eco_mammals_killerwhalesouthres'))
arcpy.Delete_management('temp')
arcpy.MakeFeatureLayer_management(ch, 'temp', """Population_EN = 'Northeast Pacific Northern Resident'""")
arcpy.CopyFeatures_management('temp', os.path.join(gdb_out, 'eco_mammals_killerwhalenorthres'))
arcpy.Delete_management('temp')



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
        os.path.join(gdb_out, 'tmp_' + rast.split('.')[0]),
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
    outname = f'eco_coarse_substrate{fc}'
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
    1
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

# NO LONGER USING THIS APPROACH
# There is now an existing mpatt version that is better.

# csv = os.path.join(root, r'01_original\harbour_seal\Harbour_Seal_Counts_from_Strait_of_Georgia_Haulout_Locations.csv')
# # some points have a positive longitude
# pts = pd.read_csv(csv)
# pts['Longitude'] = pts.Longitude.abs() * -1
# pts.to_csv(os.path.join(root, r'01_original\harbour_seal\edited.csv'))
# csv = os.path.join(root, r'01_original\harbour_seal\edited.csv')

# arcpy.env.workspace = gdb_out
# arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)
# arcpy.XYTableToPoint_management(
#     csv,
#     'temp_points',
#     'Longitude',
#     'Latitude',
#     coordinate_system = 4326
# )

# os.remove(csv)

# arcpy.Buffer_analysis('temp_points', 'temp_buffer', 200)
# arcpy.Dissolve_management(
#     'temp_buffer', 
#     'mpatt_eco_mammals_harboursealhaulout_BCMCA', # this might get merged later
#     multi_part='SINGLE_PART')

# arcpy.Delete_management('temp_points')
# arcpy.Delete_management('temp_buffer')



#######################################
#######################################
# Steller sea lion haulouts
# Rookeries and haulouts into separate datasets
# Buffer following BCMCA methodology
# 200m haulouts, 15km rookeries

# NO LONGER USING THIS CODE
# We now have a newer, already processed, mpatt dataset.

# csv = os.path.join(root, r'01_original\steller_sealion\Steller_Sea_Lion_Summer_counts_from_Haulout_Locations.csv')

# # from metadata:
# # Y: Year-round haulout site
# # W: Winter haulout site
# # A: Breeding rookeries

# # other codes in dataset: Null, ?, R, W/Y, Y/A, Y/R, Y/W
# # also, there are some records without coordinates

# df = pd.read_csv(csv)
# haul_rook = [['Y', 'W', 'W/Y', 'Y/W'], ['R', 'Y/A', 'Y/R']]
# names = ['haulout_BCMCA', 'rookery'] # haulout might be merged with another dataset
# buffs = [200, 15000]
# arcpy.env.workspace = gdb_out
# arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)

# for hr, name, buff in zip(haul_rook, names, buffs):
    
#     df_temp = df[df['SITE TYPE'].isin(hr) & ~df.LATITUDE.isnull()]
#     x = np.array(np.rec.fromrecords(df_temp.values))
#     names = df_temp.dtypes.index.tolist()
#     x.dtype.names = tuple(names)
#     arcpy.da.NumPyArrayToTable(x, os.path.join(arcpy.env.workspace, f'temp_table'))

#     arcpy.XYTableToPoint_management(
#         'temp_table',
#         'temp_points',
#         'LONGITUDE',
#         'LATITUDE',
#         coordinate_system = 4326
#     )

#     arcpy.Buffer_analysis('temp_points', 'temp_buffer', buff)
#     arcpy.Dissolve_management(
#         'temp_buffer', 
#         f'mpatt_eco_mammals_stellersealion{name}', 
#         multi_part='SINGLE_PART')

#     arcpy.Delete_management('temp_points')
#     arcpy.Delete_management('temp_table')
#     arcpy.Delete_management('temp_buffer')



#######################################
#######################################
# Fjords
# Select from oceanograhic regions dataset

# DOESN'T OVERLAP the study area. No longer needed.

# arcpy.env.workspace = gdb_out
# fc = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\mpatt_eco_coarse_oceanographicRegions_data')
# where = "Id = 11"
# arcpy.MakeFeatureLayer_management(
#     fc, 
#     'temp_lyr', 
#     where
#     )
# arcpy.CopyFeatures_management('temp_lyr', 'mpatt_eco_coarse_fjords')
# arcpy.Delete_management('temp_lyr')



#######################################
#######################################
# Kelp biobands
# buffer, dissolve, merge, clip

canopy = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\bc_ShoreZone_canopy_kelps_march2020')
under = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\bc_ShoreZone_under_kelp_march2020')
arcpy.env.workspace = gdb_out

# Arc might crash on these. Just keep rerunning.

arcpy.Buffer_analysis(canopy, 'temp_buff_can', 1, line_end_type='FLAT')
arcpy.Buffer_analysis(under, 'temp_buff_und', 1, line_end_type='FLAT')
arcpy.Dissolve_management('temp_buff_can', 'temp_dissolve_can', multi_part='SINGLE_PART')
arcpy.Dissolve_management('temp_buff_und', 'temp_dissolve_und', multi_part='SINGLE_PART')
arcpy.Clip_analysis('temp_dissolve_can', grid, 'eco_plants_kelpbiobandscanopy')
arcpy.Clip_analysis('temp_dissolve_und', grid, 'eco_plants_kelpbiobandsunderstory')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)



#######################################
#######################################
# Eelgrass polygons

eg_poly = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\eelgrass_BC_polygons_explode')
arcpy.env.workspace = gdb_out

# polygons have significant errors in them
arcpy.CopyFeatures_management(eg_poly, 'temp_copy')
arcpy.RepairGeometry_management('temp_copy')
arcpy.CheckGeometry_management('temp_copy')

arcpy.PairwiseDissolve_analysis('temp_copy', 'temp_dissolve', multi_part='SINGLE_PART') # crashed if using normal dissolve
arcpy.PairwiseClip_analysis('temp_dissolve', grid, 'eco_plants_eelgrass')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)
arcpy.Delete_management('temp_copy_CheckGeometry')



#######################################
#######################################
# Eelgrass biobands

eg_line = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\bc_ShoreZone_eelgrass_march2020')
arcpy.env.workspace = gdb_out
arcpy.Buffer_analysis(eg_line, 'temp_buff', 1, line_end_type='FLAT')
arcpy.PairwiseDissolve_analysis('temp_buff', 'temp_dissolve', multi_part='SINGLE_PART')
arcpy.PairwiseClip_analysis('temp_dissolve', grid, 'eco_plants_eelgrassbiobands')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)


#######################################
#######################################
# Surfgrass biobands
# buffer

sg_line = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\mpatt_eco_plants_surfgrass_biobands_data')
arcpy.env.workspace = gdb_out

arcpy.Buffer_analysis(sg_line, 'temp_buff', 1, line_end_type='FLAT')
arcpy.Dissolve_management('temp_buff', 'temp_dissolve', multi_part='SINGLE_PART')
arcpy.Clip_analysis('temp_dissolve', grid, 'eco_plants_surfgrassbiobands')

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

# There are memory/size constraints when going from raster to numpy.
# If I end up needing to calculate the percentile over the entire grid, then
# I will need to break it up into individual rasters and bring them into numpy:
# https://community.esri.com/t5/python-questions/how-to-calculate-the-percentile-for-each-cell-from/td-p/130653
# I just need the 80% value, so after that, I can then just apply it to the one
# study area raster.


rugs = os.path.join(dir_in, 'DST_PilotEcologicalData_BP.gdb')
arcpy.env.workspace = rugs

# combine regional rasters
from arcpy.ia import *
from arcpy.sa import *
rast1 = arcpy.Raster('rugosity_shelfsalishsea')
rast2 = arcpy.Raster('rugosity_nearshore_WCVI')
rast_reg = arcpy.ia.Merge([rast1, rast2], 'MEAN')

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

# Get just the southern shelf bioregion
rast_all = arcpy.Raster('temp_merge')
outExtractByMask = ExtractByMask(rast_all, grid)
outExtractByMask.save('temp_clip')


# ISSUE: the Salish Sea values all seem to be quite high, so once we combine
# the top 20% are disproportionaly in the Salish Sea. We can't really identify
# why, since the methodology to produce the rugosity values is the same, but
# there is clearly something wrong.
# However, to move forward, we are going to select the top 20% of values from
# the Salish Sea, and the top 20% from the remainder of the SSB. This isn't
# perfect, but oh well.
# Extract out where the Salish Sea raster overlaps the SSB
arcpy.env.workspace = rugs
outrast1 = ExtractByMask(rast1, grid)
arcpy.env.workspace = gdb_out
outrast1.save('temp_clip_salishsea')
arcpy.env.extent = grid
arcpy.CopyRaster_management('temp_clip_salishsea', 'temp_clip_ss')

# turn full raster to null where it overlaps the salish sea raster
rast_wcvi = arcpy.Raster('temp_clip')
rast_sase = arcpy.Raster('temp_clip_ss')
outisnull = IsNull(rast_sase)
outsetnull = Con(outisnull, rast_wcvi)
outsetnull.save('temp_clip_wcvi')


# So now calculate the percentile values separately:

# It looks like the only way to properly calculate percentiles is to use numpy. 
rast_clip_ss = arcpy.Raster('temp_clip_ss')
rast_clip_wcvi = arcpy.Raster('temp_clip_wcvi')
arr_ss = arcpy.RasterToNumPyArray(rast_clip_ss)
arr_wcvi = arcpy.RasterToNumPyArray(rast_clip_wcvi)
# convert 0 values to nan so that they aren't considered (I was trying to do 
# this with masking before, but the nanpercentile doesn't work properly with
# that, even though it appears that it does).
arr_nan_ss = np.where(arr_ss == 0, np.nan, arr_ss)
arr_nan_ss = np.where(arr_nan_ss > 10, np.nan, arr_nan_ss) # also an issue with big values in this one that for some reason don't get recognized as nan values in arcgis
arr_nan_wcvi = np.where(arr_wcvi == 0, np.nan, arr_wcvi)
n_80_ss = np.nanpercentile(arr_nan_ss, 80) # this is the 80th percentile value
n_80_wcvi = np.nanpercentile(arr_nan_wcvi, 80)
# now, get just the values greater than that
arr_80_ss = np.where(arr_nan_ss < n_80_ss, np.nan, arr_nan_ss)
arr_80_wcvi = np.where(arr_nan_wcvi < n_80_wcvi, np.nan, arr_nan_wcvi)
# convert nan to 0
arr_out_ss = np.nan_to_num(arr_80_ss)
arr_out_wcvi = np.nan_to_num(arr_80_wcvi)

outRas_ss = SetNull(rast_clip_ss < n_80_ss, rast_clip_ss)
outRas_wcvi = SetNull(rast_clip_wcvi < n_80_wcvi, rast_clip_wcvi)
outRas_ss.save('temp_80percentile_ss')
outRas_wcvi.save('temp_80percentile_wcvi')
outRas_ss = SetNull(rast_clip_ss < n_80_ss, 1)
outRas_wcvi = SetNull(rast_clip_wcvi < n_80_wcvi, 1)
outRas_ss.save('temp_80percentileNoValues_ss')
outRas_wcvi.save('temp_80percentileNoValues_wcvi')

rast_all = arcpy.ia.Merge([outRas_ss, outRas_wcvi], 'FIRST')
rast_all.save('temp_merge_ALL')

# convert to polys
arcpy.RasterToPolygon_conversion(
    'temp_merge_ALL',
    'eco_coarse_highrugosity',
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
    'eco_areas_coldseeps',
    multi_part='SINGLE_PART'
)
arcpy.Delete_management('memory')




#######################################
#######################################
# Geomorphic units
# split by attribute (geogreene)

geomorph_units = os.path.join(dir_in, 'DST_PilotEcologicalData_BP.gdb/BC_Geomorphic_Units')
arcpy.env.workspace = gdb_out

# clip to grid
arcpy.Clip_analysis(geomorph_units, grid, 'temp_clip')

# view list of unique values
geogreene = []
with arcpy.da.SearchCursor('temp_clip', ['GeoGreene']) as cursor:
    for row in cursor:
        if row[0] not in geogreene:
            geogreene.append(row[0])

atts = {
    'fjordwall': 'Fjord, Wall, steeply sloping',
    'fjorddepression': 'Fjord, Depression',
    'fjordcrest': 'Fjord, Crest',
    'fjordmound': 'Fjord, Mound',
    'fjorddepressionfloor': 'Fjord, Depression floor',
    'slopecanyonfloor':'Slope, Canyon floor',
    'sloperidge':'Slope, Ridge',
    'shelfwallsloping':'Shelf, Wall, sloping',
    'shelfcrest':'Shelf, Crest',
    'shelfdepressionfloor':'Shelf, Depression floor',
    'shelfmound':'Shelf, Mound',
    'shelfdepression':'Shelf, Depression',
    'shelfdepressionfloor':'Shelf, Depression floor'
}

for key in atts:

    where = f""""GeoGreene" = '{atts[key]}'"""
    arcpy.MakeFeatureLayer_management(
        'temp_clip',
        'temp_out',
        where_clause=where
    )
    out_name = f'eco_coarse_geomorphicunits_{key}'
    arcpy.CopyFeatures_management('temp_out', out_name)
    arcpy.Delete_management('temp_out')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)


#######################################
#######################################
# Marine bird colonies

gdb_birds = os.path.join(dir_in, 'mpatt_eco_birds_colonies.gdb')
arcpy.env.workspace = gdb_birds
fcs = arcpy.ListFeatureClasses()

for fc in fcs:
    
    # delete if no features intersected with the grid
    arcpy.Clip_analysis(fc, grid, 'memory/clip')
    count = arcpy.GetCount_management('memory/clip')
    if int(count[0]) == 0:
        print(f'{fc}, no features overlapping')
        arcpy.Delete_management('memory')
        continue
    bird = fc.split('_')[3]
    out_name = os.path.join(gdb_out, f'eco_birds_{bird}_colonies')
    arcpy.Dissolve_management('memory/clip', out_name, multi_part='SINGLE_PART')
    arcpy.Delete_management('memory')




#######################################
#######################################
# Important areas
# Each group requires difference geoprocessing so must be done separately

ia_dir = os.path.join(dir_in, 'ImportantAreas_Databases')
arcpy.env.workspace = gdb_out

## Birds ## 
# not including for now 


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
    if spec == 'reskillerwhl':
        continue # we are getting this one from the critical habitat dataset
    merge_list = []
    for fc in fcs:
        if spec in fc:
            merge_list.append(fc)
    arcpy.Merge_management(merge_list, 'memory/merge')
    arcpy.Dissolve_management('memory/merge', f'eco_mammals_{spec}', multi_part='SINGLE_PART')
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
arcpy.Dissolve_management('memory/merge', f'eco_inverts_spongecoral', multi_part='SINGLE_PART')
arcpy.Delete_management('memory')


## fish ##
# changed one name manually (sixgill-ed shark)
arcpy.env.workspace = os.path.join(ia_dir, 'ia_fish_wcvi.gdb')
donotinclude = ['deepseaskate', 'alaskaskate', 'deepseasole']
for fc in arcpy.ListFeatureClasses():
    species = fc.split('_')[-2]
    if species not in donotinclude:
        arcpy.CopyFeatures_management(fc, os.path.join(gdb_out, f'eco_fish_{species}'))


## invertebrates ##

# split out bivalves first
arcpy.env.workspace = os.path.join(ia_dir, 'ia_invertebrates_wcvi.gdb')
bv_wvi = 'ia_bivalves_wcvi'
arcpy.MakeFeatureLayer_management(bv_wvi, 'temp', """species = 'Razor Clam'""")
arcpy.CopyFeatures_management('temp', 'ia_razorclam_wcvi')
arcpy.Delete_management('temp')
# arcpy.MakeFeatureLayer_management(bv_wvi, 'temp', """species IN ('Olympia oyster', 'Pacific oyster')""")
# arcpy.CopyFeatures_management('temp', 'ia_oyster_wcvi')
# arcpy.Delete_management('temp')
arcpy.MakeFeatureLayer_management(bv_wvi, 'temp', """species = 'Olympia oyster'""")
arcpy.CopyFeatures_management('temp', 'ia_oysterolympia_wcvi')
arcpy.Delete_management('temp')

# split out shrimp
bv_wvi = 'ia_shrimp_wcvi'
arcpy.MakeFeatureLayer_management(bv_wvi, 'temp', """species IN ('Pandalopsis dispar, Pandalus jordani', 'Pandalus jordani')""")
arcpy.CopyFeatures_management('temp', 'ia_shrimpsmoothpink_wcvi')
arcpy.Delete_management('temp')
arcpy.MakeFeatureLayer_management(bv_wvi, 'temp', """species = 'Pandalopsis dispar, Pandalus jordani'""")
arcpy.CopyFeatures_management('temp', 'ia_shrimpsidestripe_wcvi')
arcpy.Delete_management('temp')
arcpy.MakeFeatureLayer_management(bv_wvi, 'temp', """species = 'Pandalus borealis'""")
arcpy.CopyFeatures_management('temp', 'ia_shrimpspinypink_wcvi')
arcpy.Delete_management('temp')

for fc in arcpy.ListFeatureClasses():
    species = fc.split('_')[-2]
    if species != 'bivalves' and species != 'shrimp':
        arcpy.CopyFeatures_management(fc, os.path.join(gdb_out, f'eco_inverts_{species}'))
arcpy.Delete_management('ia_razorclam_wcvi')
arcpy.Delete_management('ia_oysterolympia_wcvi')
arcpy.Delete_management('ia_shrimpsmoothpink_wcvi')
arcpy.Delete_management('ia_shrimpsidestripe_wcvi')
arcpy.Delete_management('ia_shrimpspinypink_wcvi')


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

arcpy.env.workspace = gdb_out
for spec in species:
    merge_list = []
    for fc in fcs:
        if spec in fc:
            merge_list.append(fc)
    arcpy.Merge_management(merge_list, 'memory/merge')
    if spec == 'leatherbackseaturtle':
        outname = f'eco_reptiles_{spec}'
    elif spec == 'northernfurseal':
        outname = f'eco_mammals_{spec}'
    elif spec == 'seaotter':
        outname = f'eco_mammals_{spec}_TEMP' # will get merged
    else:
        arcpy.Delete_management('memory')
        continue
    print(spec)
    print(outname)
    arcpy.Dissolve_management('memory/merge', outname, multi_part='SINGLE_PART')
    arcpy.Delete_management('memory')



#######################################
#######################################
# ADDITIONAL PROCESSING to deal with some of the above datasets that are
# duplicates and need to be merged.

# Otter datasets - merge
arcpy.env.workspace = gdb_out
otter_modeled = os.path.join(dir_in, r'DSTpilot_ecologicalData.gdb\mpatt_eco_mammals_seaotter_modeledhabitat_data')
otter_ia = 'eco_mammals_seaotter_TEMP'
arcpy.Merge_management([otter_ia, otter_modeled], 'memory/temp_otter_merge')
arcpy.Dissolve_management('memory/temp_otter_merge', 'memory/temp_otter_diss', multi_part='SINGLE_PART')
arcpy.MultipartToSinglepart_management('memory/temp_otter_diss', 'eco_mammals_seaotter')
arcpy.Delete_management('memory')
arcpy.Delete_management(otter_ia)

# Harbour seal haulouts - merge
# NO LONGER USING THIS DATASET
# There is an existing mpatt dataset to use now.
# arcpy.env.workspace = gdb_out
# seal_bcmca = 'mpatt_eco_mammals_harboursealhaulout_BCMCA'
# seal_ia = 'mpatt_eco_mammals_harbourseal_TEMP'
# arcpy.Merge_management([seal_bcmca, seal_ia], 'memory/temp_merge')
# arcpy.Dissolve_management('memory/temp_merge', 'memory/temp_diss', multi_part='SINGLE_PART')
# arcpy.MultipartToSinglepart_management('memory/temp_diss', 'mpatt_eco_mammals_harboursealhaulout')
# arcpy.Delete_management('memory')
# arcpy.Delete_management(seal_bcmca)
# arcpy.Delete_management(seal_ia)

# Steller sea lion haulouts - merge
# NO LONGER USING THIS CODE
# We now have a newer, already processed, mpatt dataset.
# arcpy.env.workspace = gdb_out
# stell_bcmca = 'mpatt_eco_mammals_stellersealionhaulout_BCMCA'
# stell_ia = 'mpatt_eco_mammals_stellersealion_TEMP'
# arcpy.Merge_management([stell_bcmca, stell_ia], 'memory/temp_merge')
# arcpy.Dissolve_management('memory/temp_merge', 'memory/temp_diss', multi_part='SINGLE_PART')
# arcpy.MultipartToSinglepart_management('memory/temp_diss', 'mpatt_eco_mammals_stellersealionhaulout')
# arcpy.Delete_management('memory')
# arcpy.Delete_management(stell_bcmca)
# arcpy.Delete_management(stell_ia)



#######################################
#######################################
# Salt marsh biobands

sm_line = os.path.join(root, r'01_original\DSTpilot_ecologicalData.gdb\bc_ShoreZone_salt_marsh_march2020')
arcpy.env.workspace = gdb_out

arcpy.Buffer_analysis(sm_line, 'temp_buff', 1, line_end_type='FLAT')
arcpy.Dissolve_management('temp_buff', 'temp_dissolve', multi_part='SINGLE_PART')
arcpy.Clip_analysis('temp_dissolve', grid, 'eco_areas_saltmarshbiobands')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)



#######################################
#######################################
# Biophysical Units
# There are a few dataset in this gdb. I'm using Biophysical_Units_L4B. 
# I look at the documentation, and it looks like this breaks down the 
# classification a bit more than the L4A dataset. However, I'm not completely 
# sure. 
# I separated out the Biophysic attribute into separate feature classes.

biounits = os.path.join(dir_in, 'PMECS_Biophysical_Geomorphic_Units.gdb/Biophysical_Units_L4B')
arcpy.env.workspace = gdb_out

# view list of unique values
biophysic = []
with arcpy.da.SearchCursor(biounits, ['Biophysic']) as cursor:
    for row in cursor:
        if row[0] not in biophysic:
            biophysic.append(row[0])

# I could do this by just changing to lowercase, but I'll keep the dictionary
# in case we change names again.
atts = {
    'trough': 'Trough',
    'shelf': 'Shelf',
    'otherbank': 'OtherBank',
    'dogfishbank': 'DogfishBank',
    'slope': 'Slope'
}

for key in atts:

    where = f""""Biophysic" = '{atts[key]}'"""
    arcpy.MakeFeatureLayer_management(
        biounits,
        'temp_out',
        where_clause=where
    )
    out_name = f'eco_coarse_biophysicalunits_{key}'
    arcpy.CopyFeatures_management('temp_out', out_name)
    arcpy.Delete_management('temp_out')



#######################################
#######################################
# Upper Ocean Sugregions

osubreg = os.path.join(dir_in, 'bcmca_eco_physical_oceanographicregions_data/PacificCoastUpperOceanSubRegionsApril2013.shp')
arcpy.env.workspace = gdb_out

# get list of ones that intersect with SSB (just do this manually)
region_id = [14,16,17,2,21,22,23,24,6]

atts = {
    'capescotttidalmixing': 14,
    'vancouverislandshelfbreak': 16,
    'vancouverislandinnershelf': 17,
    'coastalmixingregion': 2,
    'interiorgulfislands': 21,
    'harostraitandrosariopassage': 22,
    'juandefucastrait': 23,
    'juandefucaeddy': 24,
    'lowflownearshore': 6    
}

for key in atts:

    where = f""""Id" = {atts[key]}"""
    arcpy.MakeFeatureLayer_management(
        osubreg,
        'temp_out',
        where_clause=where
    )
    out_name = f'eco_coarse_oceansubreg_{key}'
    arcpy.CopyFeatures_management('temp_out', out_name)
    arcpy.Delete_management('temp_out')
    


#######################################
#######################################
# Benthic habitat classes

bhab = os.path.join(dir_in, 'Benthic_Classes/Habitats_forDST.gdb/benthic_habitats_dissolved')
arcpy.env.workspace = gdb_out

# view list of unique values
habitats = []
with arcpy.da.SearchCursor(bhab, ['Habitat']) as cursor:
    for row in cursor:
        if row[0] not in habitats:
            habitats.append(row[0])

for hab in habitats:

    where = f""""Habitat" = '{hab}'"""
    arcpy.MakeFeatureLayer_management(
        bhab,
        'temp_out',
        where_clause=where
    )
    out_hab = ''.join(hab.split()).lower()
    out_name = f'eco_coarse_benthichabitat_{out_hab}'
    arcpy.CopyFeatures_management('temp_out', out_name)
    arcpy.Delete_management('temp_out')





################################################################################
# Human activity data
################################################################################



#######################################
#######################################
# Commercial catch data
# and salmon catch data

catch_com = os.path.join(dir_in, 'CommercialFishing/hu_commercialfishing_gfsh_20211220.gdb/all_fisheries_filtered_gridded')
arcpy.env.workspace = gdb_out

# clip to grid
arcpy.Clip_analysis(catch_com, grid, 'temp_clip')

# view list of unique values
names = []
with arcpy.da.SearchCursor('temp_clip', ['name_label']) as cursor:
    for row in cursor:
        if row[0] not in names:
            names.append(row[0])

# From Carrie:
# Bottom Trawl - trawl
# Geoduck - pressure hose
# GSU - dive
# Halibut -hook and line
# Halibut and Sablefish combo - hook and line
# Lingcod - hook and line
# Midwater trawl - trawl
# Prawn and shrimp by trap - trap
# Red sea urchin - dive
# Rockfish -hook and line
# Sablefish - trap and hook and line
# Sea Cucumber - dive
# Shrimp by Trawl - trawl

atts = {
    'trawl': ['Bottom Trawl', 'Midwater Trawl', 'Shrimp by Trawl'],
    'pressurehose': ['Geoduck'],
    'hookandline': ['Halibut', 'Halibut and Sablefish combo', 'Lingcod', 'Rockfish', 'Sablefish'],
    'dive': ['Green Sea Urchin', 'Red Sea Urchin', 'Sea Cucumber'],
    'trap': ['Prawn and Shrimp by Trap']
}

# separate out
for key in atts:

    if len(atts[key]) > 1:
        where = f"""name_label IN {tuple(atts[key])}"""
    else:
        where = f""""name_label" = '{atts[key][0]}'"""
    arcpy.MakeFeatureLayer_management(
        'temp_clip',
        'temp_out',
        where_clause=where
    )
    out_name = f'temp_{key}'
    arcpy.CopyFeatures_management('temp_out', out_name)
    arcpy.Delete_management('temp_out')

# calculate total and dissolve
for key in atts:
    in_name = f'temp_{key}'
    out_name = f'temp_{key}_dissolve'
    # add field for total_kg
    # calculate from total_kg or total_kg_filtered
    arcpy.AddField_management(in_name, 'total_kg_DST', 'FLOAT')
    with arcpy.da.UpdateCursor(in_name, ['total_kg_unfiltered', 'total_kg_filtered', 'total_kg_DST']) as cursor:
        for row in cursor:
            if row[0] > 0:
                row[2] = row[0]
            else:
                row[2] = row[1]
            cursor.updateRow(row)
    arcpy.Dissolve_management(in_name, out_name, ['PU_ID'], [['total_kg_DST', 'SUM']], 'SINGLE_PART')
    arcpy.AlterField_management(out_name, 'SUM_total_kg_DST', 'total_kg_DST', 'total_kg_DST')


### Salmon data:
# combine seine and gill with trawl, and troll with hook and line
salmon_ds = os.path.join(dir_in, 'CommercialFishing/SLPSA_Salmon_Grid.gdb/SLPSA_1kmGrid_{}')
gill = salmon_ds.format('Gill')
seine = salmon_ds.format('Seine')
troll = salmon_ds.format('Troll')

# merge Gill and Seine
arcpy.Merge_management([gill, seine], 'temp_salmon_gillseine')
arcpy.Clip_analysis('temp_salmon_gillseine', grid, 'temp_salmon_gillseine_clip')
# add field and calculate
arcpy.AddField_management('temp_salmon_gillseine_clip', 'total_kg_DST', 'FLOAT')
with arcpy.da.UpdateCursor('temp_salmon_gillseine_clip', ['kg_all', 'total_kg_DST']) as cursor:
    for row in cursor:
        row[1] = row[0]
        cursor.updateRow(row)
# merge with trawl
arcpy.Merge_management(['temp_trawl_dissolve', 'temp_salmon_gillseine_clip'], 'temp_trawl_dissolvemerge')
arcpy.Dissolve_management('temp_trawl_dissolvemerge', 'temp_trawl_dissolvemergedissolve', ['PU_ID'], [['total_kg_DST', 'SUM']], 'SINGLE_PART')
arcpy.AlterField_management('temp_trawl_dissolvemergedissolve', 'SUM_total_kg_DST', 'total_kg_DST', 'total_kg_DST')
# delete and rename so naming is consistent
arcpy.Delete_management('temp_trawl_dissolve')
arcpy.Delete_management('temp_trawl_dissolvemerge')
arcpy.Rename_management('temp_trawl_dissolvemergedissolve', 'temp_trawl_dissolve')

# troll: add field and calculate, merge with hook and line, dissolve
arcpy.Clip_analysis(troll, grid, 'temp_salmon_troll_clip')
arcpy.AddField_management('temp_salmon_troll_clip', 'total_kg_DST', 'FLOAT')
with arcpy.da.UpdateCursor('temp_salmon_troll_clip', ['kg_all', 'total_kg_DST']) as cursor:
    for row in cursor:
        row[1] = row[0]
        cursor.updateRow(row)
arcpy.Merge_management(['temp_hookandline_dissolve', 'temp_salmon_troll_clip'], 'temp_hookandline_merge')
arcpy.Dissolve_management('temp_hookandline_merge', 'temp_hookandline_mergedissolve', ['PU_ID'], [['total_kg_DST', 'SUM']], 'SINGLE_PART')
arcpy.AlterField_management('temp_hookandline_mergedissolve', 'SUM_total_kg_DST', 'total_kg_DST', 'total_kg_DST')
arcpy.Delete_management('temp_hookandline_dissolve')
arcpy.Rename_management('temp_hookandline_mergedissolve', 'temp_hookandline_dissolve')

# go through gear types again and copy those with suffix _dissolve
for fc in arcpy.ListFeatureClasses('*_dissolve'):
    gear = fc.split('_')[1]
    arcpy.CopyFeatures_management(fc, f'hu_co_fishing_{gear}')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)



#######################################
#######################################
# Sport fishing

anadromous = os.path.join(dir_in, r'sport_fishing\bcmca_hu_sportfish_anadromous_data\bcmca_hu_sportfish_anadromous_data.shp')
groundfish = os.path.join(dir_in, r'sport_fishing\bcmca_hu_sportfish_groundfish_data\bcmca_hu_sportfish_groundfish_data.shp')
crab_trap = os.path.join(dir_in, r'sport_fishing\bcmca_hu_sportfish_crab_data\bcmca_hu_sportfish_crab_data.shp')
shrimp_trap = os.path.join(dir_in, r'sport_fishing\bcmca_hu_sportfish_prawnandshrimp_data\bcmca_hu_sportfish_prawnandshrimp_data.shp')
arcpy.env.workspace = gdb_out

# merge
arcpy.Merge_management([anadromous, groundfish], 'temp_hookandline_merge')
arcpy.Merge_management([crab_trap, shrimp_trap], 'temp_trap_merge')
# dissolve
arcpy.Dissolve_management('temp_hookandline_merge', 'hu_rf_fishing_hookandline', multi_part='SINGLE_PART')
arcpy.Dissolve_management('temp_trap_merge', 'hu_rf_fishing_trap', multi_part='SINGLE_PART')

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)



#######################################
#######################################
# Floating structures

fs_gdb = os.path.join(dir_in, 'Floating_Structures_PNW/floating_infrastructure.gdb')
arcpy.env.workspace = fs_gdb
fcs = arcpy.ListFeatureClasses()

arcpy.Merge_management(fcs, 'memory/merge')
arcpy.Buffer_analysis('memory/merge', 'memory/buff', 1, dissolve_option='NONE')
arcpy.env.workspace = gdb_out
arcpy.Dissolve_management('memory/buff','hu_ot_floatingstructures', multi_part='SINGLE_PART')
arcpy.Delete_management('memory')



#######################################
#######################################
# Ports and Terminals

ports = os.path.join(dir_in, 'GSR_PORTS_TERMINALS_SVW.gdb/WHSE_IMAGERY_AND_BASE_MAPS_GSR_PORTS_TERMINALS_SVW')
arcpy.env.workspace = gdb_out
arcpy.Buffer_analysis(ports, 'memory/buff', 1, dissolve_option='NONE')
arcpy.Dissolve_management('memory/buff','hu_tr_portsandterminals', multi_part='SINGLE_PART')
arcpy.Delete_management('memory') 



#######################################
#######################################
# Anchorages

# commercial
anch_co = os.path.join(dir_in, 'CanadianAnchoragesAndAnchorageAreas/ACHBRT_P.shp')
arcpy.env.workspace = gdb_out
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)
arcpy.Buffer_analysis(anch_co, 'memory/buff', 1, dissolve_option='NONE')
arcpy.Dissolve_management('memory/buff','hu_ot_anchoragescommercial', multi_part='SINGLE_PART')
arcpy.Delete_management('memory') 

# recreational
anch_re = os.path.join(dir_in, 'CanadianAnchoragesAndAnchorageAreas/ACHARE_P.shp')
arcpy.env.workspace = gdb_out
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(3005)
arcpy.Buffer_analysis(anch_re, 'memory/buff', 1, dissolve_option='NONE')
arcpy.Dissolve_management('memory/buff','hu_ot_anchoragesrecreational', multi_part='SINGLE_PART')
arcpy.Delete_management('memory') 



#######################################
#######################################
# BC cities points

cities = os.path.join(dir_in, r'MajorCities_shp\BC_MAJOR_CITIES_with_poplndata_2011_2019.shp')
arcpy.env.workspace = gdb_out
from arcpy.sa import *

# kernel density on population values
arcpy.env.extent = grid
outKD = KernelDensity(cities, 'F2019', 1000, 30000, out_cell_values='EXPECTED_COUNTS')
outKD.save('temp_kd')
# can't convert to polygon with float values. Workaround:
outMult = outKD * 100
outInt = Int(outMult)
outInt.save('temp_kd_int')

arcpy.RasterToPolygon_conversion('temp_kd_int', 'temp_poly', 'SIMPLIFY', '', 'SINGLE_OUTER_PART')
arcpy.Clip_analysis('temp_poly', grid, 'temp_poly_clip')
# add field and calculate from value * 100
arcpy.AddField_management('temp_poly_clip', 'pop_count_DST', 'FLOAT')
with arcpy.da.UpdateCursor('temp_poly_clip', ['gridcode', 'pop_count_DST']) as cursor:
    for row in cursor:
        if row[0] == 0:
            cursor.deleteRow() # delete the zero value polygon
        else:
            row[1] = row[0]/100.0
            cursor.updateRow(row)

# assign values to grid
arcpy.env.qualifiedFieldNames = False
arcpy.SpatialJoin_analysis(grid, 'temp_poly_clip', 'hu_ot_citypopulation', 
                           'JOIN_ONE_TO_MANY', 'KEEP_COMMON', 
                           match_option='HAVE_THEIR_CENTER_IN')

for field in arcpy.ListFields('hu_ot_citypopulation'):
    if not field.required and field.name != 'pop_count_DST':
        arcpy.DeleteField_management('hu_ot_citypopulation', field.name)

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)
for fc in arcpy.ListRasters('temp*'):
    arcpy.Delete_management(fc)



#######################################
#######################################
# Tenures

tenures = os.path.join(dir_in, r'TANTALIS_Crown_Tenures\TA_CROWN_TENURES_SVW.gdb\WHSE_TANTALIS_TA_CROWN_TENURES_SVW')
arcpy.env.workspace = gdb_out

# Aquaculture
where = """TENURE_PURPOSE = 'AQUACULTURE' And TENURE_SUBPURPOSE = 'SHELL FISH'"""
arcpy.MakeFeatureLayer_management(tenures, 'temp_out', where_clause=where)
arcpy.CopyFeatures_management('temp_out', 'hu_aq_aquacultureshellfish')
arcpy.Delete_management('temp_out')
where = """TENURE_PURPOSE = 'AQUACULTURE' And TENURE_SUBPURPOSE = 'FIN FISH'"""
arcpy.MakeFeatureLayer_management(tenures, 'temp_out', where_clause=where)
arcpy.CopyFeatures_management('temp_out', 'hu_aq_aquaculturefinfish')
arcpy.Delete_management('temp_out')

# Log handling
where = """TENURE_SUBPURPOSE = 'LOG HANDLING/STORAGE'"""
arcpy.MakeFeatureLayer_management(tenures, 'temp_out', where_clause=where)
arcpy.CopyFeatures_management('temp_out', 'hu_ot_loghandlingstorage')
arcpy.Delete_management('temp_out')

# Agriculture
# where = """TENURE_PURPOSE = 'AGRICULTURE'"""
# arcpy.MakeFeatureLayer_management(tenures, 'temp_out', where_clause=where)
# arcpy.CopyFeatures_management('temp_out', 'hu_ot_agriculture')
# arcpy.Delete_management('temp_out')
# DOES NOT OVERLAP WITH SSB PLANNING UNITS

# Industrial
where = """TENURE_PURPOSE = 'INDUSTRIAL' And TENURE_SUBPURPOSE NOT IN ('LOG HANDLING/STORAGE')"""
arcpy.MakeFeatureLayer_management(tenures, 'temp_out', where_clause=where)
arcpy.CopyFeatures_management('temp_out', 'hu_ot_industrial')
arcpy.Delete_management('temp_out')

# Underwater infrastructure
where = """TENURE_SUBPURPOSE IN ('ELECTRIC POWER LINE', 'SCIENCE MEASUREMENT/RESEARCH', 'SEWER/EFFLUENT LINE', 'TELECOMMUNICATION LINE', 'WATER LINE')"""
arcpy.MakeFeatureLayer_management(tenures, 'temp_out', where_clause=where)
arcpy.CopyFeatures_management('temp_out', 'hu_ot_underwaterinfrastructure')
arcpy.Delete_management('temp_out')



#######################################
#######################################
# Vessel traffic

# From Lilly (email thread): 
# Attached is the shp of the number of MMSI per grid cell (MMSI_n).
# The field All_VCou_3 is the log10(MMSI_n + 1) for all rows >0.

vtraff = os.path.join(dir_in, r'vessel_traffic\MMSI_count_1KmGrid_project.shp')
arcpy.env.workspace = gdb_out

# create count field, delete rows that are zero
arcpy.Clip_analysis(vtraff, grid, 'temp_clip')
arcpy.AddField_management('temp_clip', 'mmsi_n_DST', 'FLOAT')
with arcpy.da.UpdateCursor('temp_clip', ['All_VCou_3', 'mmsi_n_DST']) as cursor:
    for row in cursor:
        if row[0] == 0:
            cursor.deleteRow() # delete the zero value polygon
        else:
            row[1] = row[0]
            cursor.updateRow(row)

# assign values to grid
arcpy.env.qualifiedFieldNames = False
arcpy.SpatialJoin_analysis(grid, 'temp_clip', 'hu_tr_vesseltraffic', 
                           'JOIN_ONE_TO_MANY', 'KEEP_COMMON', 
                           match_option='HAVE_THEIR_CENTER_IN')

for field in arcpy.ListFields('hu_tr_vesseltraffic'):
    if not field.required and field.name != 'mmsi_n_DST':
        arcpy.DeleteField_management('hu_tr_vesseltraffic', field.name)

for fc in arcpy.ListFeatureClasses('temp*'):
    arcpy.Delete_management(fc)


