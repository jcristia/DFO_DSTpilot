# create feature classes of Marxan results

import arcpy
import os
import pandas as pd
import numpy as np


root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot'
output_dir = 'scripts/dst_03_marxan'
gdb = 'spatial/04_outputs/marxan.gdb'
pus = os.path.join(root, r'spatial\03_working\dst_grid.gdb\dst_grid1km_PUVSP')
arcpy.env.workspace = os.path.join(root, gdb)
arcpy.env.qualifiedFieldNames = False


dirs = os.listdir(os.path.join(root, output_dir))
dirs = [dir for dir in dirs if 'scenario' in dir] # scenario needs to be in the name

for dir in dirs:

    if arcpy.Exists(dir):
        continue

    output_folder = os.path.join(root, output_dir, dir, 'output')
    best_sol = pd.read_csv(os.path.join(output_folder, 'output_best.csv'))

    # copy pus to gdb
    arcpy.CopyFeatures_management(pus, dir)

    # add field
    arcpy.AddField_management(dir, 'solution_zone', 'SHORT')

    # with update cursor, for planning unit get zone id minus 1
    with arcpy.da.UpdateCursor(dir, ['uID', 'solution_zone']) as cursor:
        for row in cursor:
            sz = best_sol.zone[best_sol.planning_unit == row[0]] - 1
            row[1] = sz.values[0]
            cursor.updateRow(row)


