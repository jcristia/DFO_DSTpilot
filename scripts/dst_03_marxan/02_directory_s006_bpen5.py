# create Marxan directory and input.dat for a new scenario

import shutil
import pandas as pd
import arcpy


#####################################
# General inputs
# see each individual file section for additional inputs to potentially change
name = 'scenario_006_bpen5'

# set mpas to zero cost
# For lockedin scenarios, I am not changing the cost of MPA pus to zero. I want 
# to see how much more efficient a solution is when we don't have to lock in 
# mpas. I am only setting them to zero when I have boundary penalties (which, I
# guess will be every solution when working in Marxan).
mpas_to_zero = True

to_not_include = [ # optional files to not include in input.dat
    'puzone', # no set up for this yet
    'zonetarget2', # no set up for this yet
    'zonecontrib', # no set up for this yet
    'zonecontrib2', # no set up for this yet
    ]

# For our scenarios we are only using: PULOCK, BOUND, ZONETARGET, ZONEBOUNDCOST  
# OPTIONAL FILES:
# see the manual starting on page 5-58 for descriptions
# 'bound' # needed if BLM is set to 1
# 'pulock', # for locking in pus to just 1 zone each (can be different zones)
# 'puzone', # for locking in pus to multiple zones (i.e. can be assigned to one zone or another). Do not use if pulock is being used.
# 'zonetarget', # set targets by zone
# 'zonetarget2', # simplified alternative to zonetarget
# 'zoneboundcost', # Controls the BLM within a zone and also to prevent/enable certain zones to border each other
# 'zonecontrib', # sets up proportion that a pu in a zone can contribute to meeting targets (not related to target values though)
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
# create directory
if not os.path.isdir(scenario):
    os.mkdir(scenario)
    os.mkdir(os.path.join(root, scenario, 'input'))
    os.mkdir(os.path.join(root, scenario, 'output'))

# copy in executable and input.dat file
executable = os.path.join(marxan_exe, 'MarZone_x64.exe')
shutil.copy2(executable, scenario)


#####################################
# copy in and configure input.dat file
# I'm only controling the BLM value here. For any other values, just edit it
# manually.

BLM = 1 #set to 1 or 0. Actual BLM is set in zoneboundcost if this is 1.

inputdat_orig = os.path.join(marxan_exe, 'input_TEMPLATE.dat')
inputdat = os.path.join(root, scenario, 'input_TEMPLATE.dat')
shutil.copy(inputdat_orig, inputdat)

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

if mpas_to_zero:
    field_names = ['uID']
    cursor = arcpy.da.SearchCursor(mpas, field_names)
    df_mpas = pd.DataFrame(data=[row for row in cursor], columns=field_names)
    df = df.merge(df_mpas, how='left', left_on='id', right_on='uID')
    df['cost'] = df['cost'].where(df.uID.isna(), 0)
    df = df[['id', 'cost']]

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
# Create zone targets file (zonetarget.dat)
# This is set up assuming that feature targets are all met in only one zone. If
# that ever changes then we will need to update this.

if not 'zonetarget' in to_not_include:

    # zoning csv
    zone_column = 'zone_20220315_4zones'
    zones = pd.read_csv(zoning)
    # feature ids
    features = pd.read_csv(os.path.join(root, scenario, 'input/feat.dat'), sep='\t')
    # marxan zone ids
    zoneids = pd.read_csv(os.path.join(root, scenario, 'input/zones.dat'), sep='\t')

    df_all = pd.DataFrame(columns=['zoneid', 'featureid', 'target', 'targettype'])
    for z in zoneids.zoneid[1:]: # we skip the first 'available' Marxan zone
        
        # the original zone id in the zones cvs
        z_original = z - 1
        
        # the feature names in each zone
        fzones = zones.feature[zones[zone_column]==z_original]

        # with those names, get the feature ids and targets
        fids = features[features.name.isin(fzones)]
        fids = fids[['id', 'prop']]
        fids = fids.rename(columns={'id':'featureid', 'prop':'target'})

        fids['zoneid'] = z # marxan zone id
        fids['targettype'] = 1
        fids = fids[['zoneid', 'featureid', 'target', 'targettype']]
        df_all = df_all.append(fids)

    # This file for some reason has to be sorted from lowest to highest value, first
    # by zoneid, then by featureid, and then by target.
    df_all = df_all.sort_values(by=['zoneid', 'featureid'])

    df_all.to_csv(os.path.join(root, scenario, 'input/zonetarget.dat'), index=False, sep='\t')


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
    'cost':   [0,0,0,0,0,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001]
    }
    df_zbc = pd.DataFrame.from_dict(zbc)
    df_zbc.to_csv(os.path.join(root, scenario, 'input/zoneboundcost.dat'), index=False, sep='\t')




#####################################
# Create boundary length file (bound.dat)

# This file lists the boundary costs between 2 planning units. In our case, 
# these will all be 10 since we are using modified pu "areas" of 100.

# I am modifying the code from the ArcMarxan toolbox. I want a standalone 
# script. I don't want to rely on a GUI in Arcmap 10.x.
# source: https://www.aproposinfosystems.com/en/solutions/arcgis-plugins/arcmarxan-toolbox/

# In general, when thinking about boundary "costs":
# Imagine selecting a bunch of units where none are connected. The total cost is
# just the cost of the units. I think the boundary penalty stuff works by 
# subtracting away from this the shared boundaries so that it rewards solutions
# that are more clustered. Therefore, our total cost will still never be more 
# than just the total cost of the units.

# Why do we give a shared boundary value to a cell matched to itself?
# (I think...) an effecient solution is very clumped and therefore a lot of
# boundary gets subtracted away to lower the total cost. However, normally a 
# cell on the edge will have 1-3 boundaries that aren't shared with another 
# cell, so there is nothing to subtract away if it gets selected. So, we add 
# self connections to itself to artificaially give it something to subtract away
# so that it has a chance to get selected in a solution that is prioritizing a
# minimized boundary. However, we don't want to make it too desirable, so we 
# give it just half the amount(??).

def boundary(boundary_layer, pu_field, boundary_treatment, weighting):

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
    temp_gdb = os.path.join(scenario, 'temp.gdb')
    
    # buffer and dissolve layer in one step
    buffFeatClass = os.path.join(temp_gdb, "amt_temp_buffer")
    arcpy.Buffer_analysis(boundary_layer, buffFeatClass, '100 Meters', dissolve_option="ALL")
    # union buffered layer with the pulayer and delete temp source
    unionFeatClass = os.path.join(temp_gdb, "amt_temp_union")
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
    arcpy.PolygonNeighbors_analysis(unionFeatClass,outDBF,[pu_field],both_sides="NO_BOTH_SIDES",out_linear_units="METERS")
    boundList = []
    x = 0
    for row in arcpy.da.SearchCursor(outDBF,'*'):
        x += 1
        id1 = int(row[1])
        id2 = int(row[2])
        if id1 == -1:
            id1 = id2
        nodes = row[4]
        sideLen = row[3]
        # note: node > 0 means a touching corner, not a touching side.
        # JC: therefore, we are looking for relationships where there is a
        # shared length and not a touching corner to the boundary. All shared
        # boundaries will have 0 node corners touching.
        if nodes == 0 and sideLen > 0:
            bValue = sideLen * weighting
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
    oFileName = os.path.join(scenario,'input/bound.dat')
    oFile = open(oFileName,'w')
    oFile.write('id1\tid2\tboundary\n')
    for rec in boundList:
        oFile.write('%d\t%d\t%f\n' % (rec[0],rec[1],rec[2]))   
    oFile.close()
    arcpy.Delete_management(outDBF)
    arcpy.Delete_management(unionFeatClass)
    # delete temp gdb
    arcpy.Delete_management(temp_gdb)
    return


if not 'bound' in to_not_include:

    boundary(
        boundary_layer = puvsp_fc,
        pu_field = 'uID',
        boundary_treatment = 'Half Value',
        weighting = 0.01 # to get edges of 10
    )


#####################################
# Create zone contribution file (zonecontrib.dat)
# NO LONGER USING THIS, BUT MAKE SOME NOTES FOR UNDERSTANDING
# I misunderstood zone contributions. There explanation is a bit confusing. It 
# is for the level that a pu can contribute to meeting a target. For instance, 
# for an eco feature, if it is in a conservation zone then 100% of the area of a
# pu can contribute towards the target. However, perhaps you'll also consider 
# that some protection can be met in an industrial zone (e.g. since a fish may 
# not be fished there, but the activity could indirectly affect it). In this 
# case, perhaps only 50% of a pu can contribute towards meeting the target for 
# that fish. Therefore, by not specifying zonecontrib.dat, I just say that all 
# are treated equal, and since I don't target features in more than 1 zone then 
# I don't need to worry about this.