# create Marxan directory and input.dat for a new scenario

import shutil
import pandas as pd
import arcpy


#####################################
# inputs
name = 'scenario_001'
BLM = 1
to_not_include = [
    'zonetarget', 
    'zoneboundcost', 
    'zonecontrib2']  # optional files to not include in input.dat
# paths
root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot'
marxan_exe = r'C:\Users\jcristia\Desktop\MarxanWZones_v406'
feat_targs = os.path.join(root, r'scripts\dst_01_preprocess\DSTpilot_spatialData - naming_scheme.csv')
puvsp_fc = os.path.join(root, r'spatial\03_working\dst_grid.gdb\dst_grid1km_PUVSP')


#####################################
# create directory and input.dat file
scenario = os.path.join(root, 'scripts/dst_03_marxan', name)
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
field_names = ['uID', 'COST']
cursor = arcpy.da.SearchCursor(puvsp_fc, field_names)
df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df = df.rename(columns={'uID':'id', 'COST':'cost'})
out = os.path.join(root, scenario, 'input/pu.dat')
df.to_csv(out, index=False, sep='\t')


#####################################
# Create planning unit versus feature file (puvfeat.dat)

# puvsp amounts
field_names = [i.name for i in arcpy.ListFields(puvsp_fc) if i.type != 'OID']
cursor = arcpy.da.SearchCursor(puvsp_fc, field_names)
df = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df = df.rename(columns={'uID':'id', 'COST':'cost'})
out = os.path.join(root, scenario, 'input/pu.dat')
df.to_csv(out, index=False, sep='\t')

# feature ids
features = pd.read_csv(os.path.join(root, scenario, 'input/feat.dat'), sep='\t')

# set up blank df with 3 fields
# for each feature in features
# from puvsp get uIDs where greater than 0. 
# Add field for feature id.
# change field names
# Append to larger df

# output to dat