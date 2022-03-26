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
# zone, as opposed to just 1 large cluster. Then if you set the value between
# zones to something between 0-1 then you create clumping within zones. A
# negative value would prevent zones from bordering each other.

if not 'zoneboundcost' in to_not_include:
    pass
    # set up a file with all combinations of all zones and set them to 0 to start



#####################################
# Create boundary length file (bound.dat)
# Lists the boundary costs between 2 planning units. In our case, these will all
# be 10 since we are using modified pu "areas" of 100.

# This one is a bit more involved.
# For each pu do a spatial select? Do I list each combination and the opposite?

# Check the good practices handbook
# What to do about those on the border of the study area? Is this where I alter
# their cost.
# Why in their example to they include a pu boundary with itself?
# For each pu do a spatial select? Do I list each combination and the opposite? (e.g. 1-2 and 2-1)

if not 'bound' in to_not_include:
    pass