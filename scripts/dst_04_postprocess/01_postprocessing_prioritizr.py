# read prioritizr outputs into fgdb, calculate selection frequency

# OLD code had a lot of stuff on replacement importance, which I don't think
# I'll be calculating here. I'm going to comment it out for now.

import arcpy
import os
import pandas as pd
import numpy as np


root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot'
output_dir = 'scripts/dst_02_prioritizr/outputs'
gdb = 'spatial/04_outputs/prioritizr.gdb'
pus = os.path.join(root, r'spatial\03_working\dst_grid.gdb\dst_grid1km_PUVSP')
arcpy.env.workspace = os.path.join(root, gdb)
arcpy.env.qualifiedFieldNames = False


dirs = os.listdir(os.path.join(root, output_dir))

# create a feature class for each solution
for dir in dirs:
    
    if arcpy.Exists(dir):
        continue
    
    solution = os.path.join(root, output_dir, dir, 'solution.csv')
    df = pd.read_csv(solution)

    # solution columns
    cols = ['uID']
    sol_cols = [col for col in df.columns if col.startswith('solution')]
    cols.extend(sol_cols)
    df = df[cols]

    # calculate selection frequency
    # I'm not sure I need this anymore, and with zones it becomes more complicated.
    # However, maybe I will calculate multiple solutions for the final solution,
    # but then I can just add them manually.
    # df['sel_freq'] = df.loc[:, sol_cols].sum(axis=1) / len(sol_cols)
    # keep just the top solution
    # df = df[['uID', 'solution_1', 'sel_freq']]

    # If its a ZONE problem, then summarize into 1 column.
    # This requires that 'zone' is in the name of the first solution.
    # Keep in mind, this is only for if there is 1 zone problem per fc and not
    # a portfolio of multiple solutions.
    if 'zone' in sol_cols[0]:

        # this is the vectorized numpy way, but the pandas apply/lambda way is
        # a bit easier to understand.
        conditions = []
        outputs = []
        i=1
        for sol_col in sol_cols:
            condition = df[sol_col]==1
            conditions.append(condition)
            outputs.append(i)
            i+=1
        df['solution_1_summary'] = pd.Series(np.select(conditions, outputs, 0))

    # output to gdb table
    x = np.array(np.rec.fromrecords(df.values))
    names = df.dtypes.index.tolist()
    x.dtype.names = tuple(names)
    arcpy.da.NumPyArrayToTable(x, os.path.join(arcpy.env.workspace, 'tbl_TEMP'))
    # join to copy of point feature class
    jointbl = arcpy.AddJoin_management(pus, 'uID', 'tbl_TEMP', 'uID')
    sol_num = dir
    arcpy.CopyFeatures_management(jointbl, f'{sol_num}')
    arcpy.DeleteField_management(f'{sol_num}', ['OBJECTID_1', 'uID_1'])
    arcpy.Delete_management(jointbl)
    arcpy.Delete_management('tbl_TEMP')



### Irreplaceability ###

# I used to add these values to the df in the functions above. However, I don't
# think I will be calculating it except for the last solution. It is very slow.
# I now output it in a separate csv file. When I get to calculating it, I can
# add some functionality here to append it to the feature classes.