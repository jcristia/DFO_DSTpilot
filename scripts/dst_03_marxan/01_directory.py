# create Marxan directory and input.dat for a new scenario

import shutil
import pandas as pd
import arcpy


#####################################
# inputs
# see each individual file section for additional inputs to potentially change
name = 'scenario_001'
BLM = 1 #set to 1 or 0. BLM is set in zoneboundcost
to_not_include = [ # optional files to not include in input.dat
    'puzone', # no set up for this yet
    'zonetarget', # no set up for this yet
    'zonetarget2', # no set up for this yet
    'zonecontrib2', # no set up for this yet
    ]

# For our scenarios we are only using: PULOCK, BOUND, ZONECONTRIB, ZONEBOUNDCOST  
# OPTIONAL FILES:
# see the manual starting on page 5-58 for descriptions
# 'bound' # needed if BLM is set to 1
# 'pulock', # for locking in pus to just 1 zone each (can be different zones)
# 'puzone', # for locking in pus to multiple zones (i.e. can be assigned to one zone or another). Do not use if pulock is being used.
# 'zonetarget', # only use if overall targets are not set in feat.dat
# 'zonetarget2', # simplified alternative to zonetarget
# 'zoneboundcost', # Controls the BLM within a zone and also to prevent/enable certain zones to border each other
# 'zonecontrib', # sets up targets by zones as proportion of the overall target in feat.dat
# 'zonecontrib2' # Simplified version of zonecontrib.

# paths
root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot'
marxan_exe = os.path.join(root, r'scripts\dst_03_marxan\MarxanWZones_v406')
feat_targs = os.path.join(root, r'scripts\dst_01_preprocess\DSTpilot_spatialData - naming_scheme.csv')
puvsp_fc = os.path.join(root, r'spatial\03_working\dst_grid.gdb\dst_grid1km_PUVSP')
mpas = os.path.join(root, r'spatial\03_working\dst_mpas.gdb\mpas_rcas_marxan')
zoning = os.path.join(root, r'scripts\dst_02_prioritizr\Zoning frameworks - zones_features.csv')
scenario = os.path.join(root, 'scripts/dst_03_marxan', name)

#####################################
# create directory and input.dat file
if not os.path.isdir(scenario):
    os.mkdir(scenario)
    os.mkdir(os.path.join(root, scenario, 'input'))
    os.mkdir(os.path.join(root, scenario, 'output'))

# copy in executable and input.dat file
executable = os.path.join(marxan_exe, 'MarZone_x64.exe')
shutil.copy2(executable, scenario)
inputdat_orig = os.path.join(marxan_exe, 'input_TEMPLATE.dat')
inputdat = os.path.join(root, scenario, 'input_TEMPLATE.dat')
shutil.copy(inputdat_orig, inputdat)

# configure input.dat
new_input = os.path.join(root, scenario, 'input.dat')
with open(inputdat) as oldfile, open(new_input, 'w') as newfile:
    for line in oldfile:
        if BLM != 1 and 'BLM ' in line:
            newfile.write(f"BLM {str(BLM)}\n")
            continue        
        if not any(file in line for file in to_not_include):
            newfile.write(line)
# Delete template
os.remove(inputdat)


#####################################
# Create features file (feat.dat)
# This is the list of the features (eco/hu), their numerical IDs, and overall
# targets.

features = pd.read_csv(feat_targs)
spec = os.path.join(root, scenario, 'input/feat.dat')
spec_df = features[['processed_file', 'Coarse-filter feature; Low Target Range (10%)']]
spec_df['id'] = spec_df.index + 1
spec_df =  spec_df.rename(columns={'processed_file':'name', 'Coarse-filter feature; Low Target Range (10%)':'prop'})
spec_df = spec_df[['id', 'prop', 'name']]
spec_df['prop'] = round(spec_df.prop / 100, 3)
spec_df.to_csv(spec, index=False, sep='\t')


#####################################
# Create planning unit file (pu.dat)
# This is the list of the planning unit IDs and the cost of selecting each unit.
# For us, the cost is the same and is relative to area.

field_names = ['uID', 'COST']
cursor = arcpy.da.SearchCursor(puvsp_fc, field_names)
df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df = df.rename(columns={'uID':'id', 'COST':'cost'})
out = os.path.join(root, scenario, 'input/pu.dat')
df.to_csv(out, index=False, sep='\t')


#####################################
# Create planning unit versus feature file (puvfeat.dat)
# This relates the features to the planning units and gives the amount of each
# feature in the pu.

# puvsp amounts
field_names = [i.name for i in arcpy.ListFields(puvsp_fc) if i.type != 'OID']
cursor = arcpy.da.SearchCursor(puvsp_fc, field_names)
df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df = df.rename(columns={'uID':'pu'})

# feature ids
features = pd.read_csv(os.path.join(root, scenario, 'input/feat.dat'), sep='\t')

# pull and combine
df_all = pd.DataFrame(columns=['featureid', 'pu', 'amount'])
for ind in features.index:
    feat_id = features.id[ind]
    feat = features.name[ind]
    df_featpus = df[df[feat] > 0]
    df_featpus = df_featpus[['pu', feat]]
    df_featpus['featureid'] = feat_id
    df_featpus = df_featpus.rename(columns={feat:'amount'})
    df_featpus = df_featpus[['featureid', 'pu', 'amount']]
    df_all = df_all.append(df_featpus)
df_all.to_csv(os.path.join(root, scenario, 'input/puvfeat.dat'), index=False, sep='\t')


#####################################
# Create zones file (zones.dat)
# A list of the zones and their numerical IDs.

# In Marxan, zone 1, must be 'available'.
zones = {
    'zoneid':[1,2,3,4,5],
    'zonename':['available', 'protection', 'industrial_nearshore', 'industrial_offshore', 'fishing']
}
df_z = pd.DataFrame.from_dict(zones)
df_z.to_csv(os.path.join(root, scenario, 'input/zones.dat'), index=False, sep='\t')


#####################################
# Create costs file (costs.dat)
# The cost name and cost id. We will only have 1 cost type across zones.

costs = {
    'costid':[1],
    'costname':['cost']
}
df_c = pd.DataFrame.from_dict(costs)
df_c.to_csv(os.path.join(root, scenario, 'input/costs.dat'), index=False, sep='\t')


#####################################
# Create zone cost file (zonecost.dat)
# Relates the cost type to each zone. We will only have 1 cost type across zones.

zcosts = {
    'zoneid':[2,3,4,5],
    'costid':[1,1,1,1],
    'multiplier':[1,1,1,1]
}
df_zc = pd.DataFrame.from_dict(zcosts)
df_zc.to_csv(os.path.join(root, scenario, 'input/zonecost.dat'), index=False, sep='\t')


#####################################
# Create pulock file (pulock.dat)
# Lock in MPAs exclusively to the protection zone.

if not 'pulock' in to_not_include:

    field_names = ['uID']
    cursor = arcpy.da.SearchCursor(mpas, field_names)
    df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    df = df.rename(columns={'uID':'puid'})
    df['zoneid'] = 2
    df.to_csv(os.path.join(root, scenario, 'input/pulock.dat'), index=False, sep='\t')


#####################################
# Create zone contribution file (zonecontrib.dat)
# This works with the feature targets in feat.dat. It assigns a feature to a
# a zone and the fraction of the target that can be met in that zone.
# For a species, the fraction can add up to greater than 1. This allows a target
# to be met in one zone OR another, or a combination of them.
# However, in our case, a feature is only assigned to one zone and therefore we
# will only list a feature once and its fraction will be 1.

if not 'zonecontrib' in to_not_include:

    # zoning csv
    zone_column = 'zone_20220315_4zones'
    zones = pd.read_csv(zoning)
    # feature ids
    features = pd.read_csv(os.path.join(root, scenario, 'input/feat.dat'), sep='\t')
    # marxan zone ids
    zoneids = pd.read_csv(os.path.join(root, scenario, 'input/zones.dat'), sep='\t')

    df_all = pd.DataFrame(columns=['zoneid', 'specid', 'fraction'])
    for z in zoneids.zoneid[1:]: # we skip the first 'available' Marxan zone
        
        # the original zone id in the zones cvs
        z_original = z - 1
        
        # the feature names in each zone
        fzones = zones.feature[zones[zone_column]==z_original]

        # with those names, get the feature ids
        fids = features[features.name.isin(fzones)]
        fids = fids[['id']]
        fids = fids.rename(columns={'id':'specid'})

        fids['zoneid'] = z # marxan zone id
        fids['fraction'] = 1 # all features should be represented in only 1 zone
        fids = fids[['zoneid', 'specid', 'fraction']]
        df_all = df_all.append(fids)

    df_all.to_csv(os.path.join(root, scenario, 'input/zonecontrib.dat'), index=False, sep='\t')


#####################################
# Create zone boundary cost file (zoneboundcost.dat)

# see page 5-64 in the manual details and for calibrating.
# In brief: a cost of zero means no clustering is considered. In their examples
# they use the zone to zone cost to create clumping WITHIN a zone. This is
# different than I've thought of it in the past.
# You can probably control it in different ways, but if you keep the within zone
# clumping set to zero, then there can be different clumps across space for that
# zone, as opposed to just 1 large cluster per zone. Then if you set the value 
# between zones to something between 0-1 then you create clumping within zones. 
# A negative value would prevent zones from bordering each other.
# Another way to think about it:
# By setting a value greater than zero between different zones, you are saying
# that you dont want pus of one zone embedded in the cluster of another zone.

if not 'zoneboundcost' in to_not_include:
    
    zbc = {
    'zoneid1':[1,2,3,4,5,1,1,1,1,2,2,2,3,3,4],
    'zoneid2':[1,2,3,4,5,2,3,4,5,3,4,5,4,5,5],
    'cost':   [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    }
    df_zbc = pd.DataFrame.from_dict(zbc)
    df_zbc.to_csv(os.path.join(root, scenario, 'input/zoneboundcost.dat'), index=False, sep='\t')




#####################################
# Create boundary length file (bound.dat)

# I am using the code from the ArcMarxan toolbox. I want a standalone script.
# I don't want to rely on a GUI in Arcmap 10.x.
# source: https://www.aproposinfosystems.com/en/solutions/arcgis-plugins/arcmarxan-toolbox/

# Lists the boundary costs between 2 planning units. In our case, these will all
# be 10 since we are using modified pu "areas" of 100.

# TO CLEAN UP:
'''
	- Understand why they give a boundary value to a cell matched to itself for just some cells. (the best practices manual has more on this. 
		- we add in a self connection because normally when a unit gets selected in a solution, we add up the shared boundaries and subtract it from the total shared boundaries. An effecient soluton is very clumped and therefore a lot of boundary gets subtracted away (it gets added to the total cost, and we are trying to minimize the cost). However, normally a cell on the edge will have 1-3 boundaries that aren't shared with another cell, so there is nothing to subtract away if it gets selected. So then we add self connections to itself to artificaially give it something to subtract away so that it has a chance to get selected. However, we don't want to make it too desirable, so we just give it half the amount.
		- so therefore, in the bound.dat file, there should be some self connections that are also 1000 and 1500. See if I can find these.
		- Also, I need to scale these to be 10 or 5.
		- Also, when thinking of how this works with cost. Imagine selecting a bunch of units where none are connected. The total cost is just the cost of the units. The penalty stuff works by subtracting away from this the shared boundaries. So its not like we start by adding the total possible boundary. Our total cost will still never be more than just the total cost of the units.
'''

boundary_layer = puvsp_fc
pu_field = 'uID'
boundary_method = 'Measured'
boundary_treatment = 'Half Value'
weighting = 0.1

def boundary(boundary_layer):

    #
    # The approach is straightforward. Marxan must account for outer boundaries therefore
    # a new layer must be created which has an outer boundary. The natural way to do this is
    # to dissolve the pu layer, then buffer it and union the result with the original pu layer.
    # The next step is to set the PUID of the outer buffer to -1 as -1 is not a permitted pu id 
    # in Marxan. Then the arcpy polygon neighbors analysis runs and gives us a table. When one of 
    # the puid's equals -1 then that indicates an outer boundary and the id is set to match the 
    # other id and thus becomes the "self" boundary. The end result is sorted and written to disk.
    #

    # create temp gdb for workspace
    arcpy.CreateFileGDB_management(scenario, 'temp.gdb')
    arcpy.env.workspace = os.path.join(scenario, 'temp.gdb')
    
    # buffer and dissolve layer in one step
    buffFeatClass = "amt_temp_buffer"
    arcpy.Buffer_analysis(boundary_layer, buffFeatClass, '100 Meters', dissolve_option="ALL")
    # union buffered layer with the pulayer and delete temp source
    unionFeatClass = "amt_temp_union"
    arcpy.Union_analysis([buffFeatClass, boundary_layer], unionFeatClass)
    arcpy.Delete_management(buffFeatClass)
    # update planning unit id to equal -1 
    # based on knowledge that FID_boundary_layer field will be set to -1 or 0 because 
    # the field didn't exist in the source layer for the boundary layer
    field_names = [f.name for f in arcpy.ListFields(unionFeatClass)]
    with arcpy.da.UpdateCursor(unionFeatClass,[field_names[5],pu_field]) as cursor:
        for row in cursor:
            if row[0] == -1 or row[0] == 0:
                row[1] = -1
                cursor.updateRow(row)
                break

    # get boundaries and lengths
    outDBF = os.path.join(scenario,'tempBound.dbf')
    if boundary_method in ["Weighted","Field"]:
        arcpy.PolygonNeighbors_analysis(unionFeatClass,outDBF,[pu_field,calc_field_name],both_sides="NO_BOTH_SIDES",out_linear_units="METERS")
    else:
        arcpy.PolygonNeighbors_analysis(unionFeatClass,outDBF,[pu_field],both_sides="NO_BOTH_SIDES",out_linear_units="METERS")
    boundList = []
    x = 0
    for row in arcpy.da.SearchCursor(outDBF,'*'):
        x += 1
        id1 = int(row[1])
        id2 = int(row[2])
        if id1 == -1:
            id1 = id2
        if boundary_method in ["Weighted","Field"]:
            nodes = row[6]
            sideLen = row[5]
        else:
            nodes = row[4]
            sideLen = row[3]
        # note: node > 0 means a touching corner, not a touching side
        if nodes == 0 and sideLen > 0:
            if boundary_method in ["Weighted","Field"]:
                if id1 == id2 and int(row[1]) == -1:
                    # note this is for the perimeter where the buffered calc field value is zero 
                    # and so the difference method will fail
                    fValue = row[4]
                else:
                    if field_calc_method == 'Mean':
                        fValue = (row[3]+row[4])/2.0
                    elif field_calc_method == 'Maximum':
                        if row[3] > row[4]:
                            fValue = row[3]
                        else:
                            fValue = row[4]
                    else:
                        if row[3] < row[4]:
                            fValue = row[3]
                        else:
                            fValue = row[4]
            if boundary_method == 'Measured':
                bValue = sideLen
            elif boundary_method == 'Weighted':
                bValue = sideLen * fValue
            elif boundary_method == 'Field':
                bValue = fValue
            elif boundary_method == 'Single Value':
                bValue = boundary_value
            # deal with boundary units special cases
            if id1 == id2:
                if boundary_treatment in ["Full Value","Half Value"]:
                    if boundary_treatment == "Half Value":
                        bValue = bValue / 2.0
                    boundList.append([int(id1),int(id2),bValue])
            else:
                boundList.append([int(id1),int(id2),bValue])
    boundList.sort()

    # convert to bound.dat file
    oFileName = os.path.join(marxan_input_folder,'bound.dat')
    oFile = open(oFileName,'w')
    oFile.write('id1\tid2\tboundary\n')
    for rec in boundList:
        oFile.write('%d\t%d\t%f\n' % (rec[0],rec[1],rec[2]))   
    oFile.close()
    arcpy.Delete_management(outDBF)
    arcpy.Delete_management(unionFeatClass)


    # delete temp gdb
    # FUCK ME, BE CAREFUL

    return

if not 'bound' in to_not_include:
    pass