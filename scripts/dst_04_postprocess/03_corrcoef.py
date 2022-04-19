# Calculate the Pearson correlation coefficient between two rasters.
# We want to see how similar two solutions are.

import arcpy
import numpy as np
import pandas as pd
import seaborn as sns


root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot'

# Prioritizr solution
prsol = os.path.join(root, r'spatial\04_outputs\prioritizr.gdb\s102_minset_adjtarg')

# Marxan solutions
sol_range = ['101', '102', '103', '104', '105', '106']
iterations = [ #getting a little lazy and just hardcoding this
    100000,
    1000000,
    10000000,
    100000000,
    1000000000,
    2000000000
]
marxan_gdb = os.path.join(root, r'spatial\04_outputs\marxan.gdb')
arcpy.env.workspace = marxan_gdb

# To raster to numpy matrix - prioritizr
arcpy.PolygonToRaster_conversion(
    prsol,
    'solution_1_summary',
    'memory/temp_pr',
    cellsize=1000
)
arr_pr = arcpy.RasterToNumPyArray('memory/temp_pr', nodata_to_value=9)
arr_pr = np.ma.masked_where(arr_pr==9, arr_pr)
# have to mask, or else it finds correlation between the nodata values
arr_pr = arr_pr.flatten() 
# it can only be 1 dimension
# the np.corrcoef does a row-wise comparison
# if I include them in square brackets then they can't exceed 1-d.
# outside of brackets they can be 2 dimensions but then it compares row to row


arrays = [arr_pr]

# To raster to numpy matrix - marxan
for sol in sol_range:
    mxsol = arcpy.ListFeatureClasses(f'scenario_{sol}*')[0]
    arcpy.PolygonToRaster_conversion(
    mxsol,
    'solution_zone',
    'memory/temp_mx',
    cellsize=1000
    )
    arr_mx = arcpy.RasterToNumPyArray('memory/temp_mx', nodata_to_value=9)
    arr_mx = np.ma.masked_where(arr_mx==9, arr_mx)
    arr_mx = arr_mx.flatten()
    arrays.append(arr_mx)
    arcpy.Delete_management('memory/temp_mx')

# calc corrcoef
cc_m = np.ma.corrcoef(arrays)
# I've confirmed that this approach get the same result as calculating it with
# the Band Statistics tool in Arc.

# create dataframe: number of iterations and corrcoef
ccs = list(cc_m[1:,0])
df = pd.DataFrame(list(zip(iterations, ccs)), columns=['iteration', 'corrcoef'])

# plot corrcoef vs iterations
df = df.astype({'iteration':'int64', 'corrcoef':'float32'})
sns.set()
sns.set_style('white')
sns.set_context('notebook', font_scale=1.25, rc={"lines.linewidth": 2})
g = sns.lineplot(data=df, x='iteration', y='corrcoef', marker='o', markersize=8)
#g.set(ylim=(0.91, 1.01))
g.set(xscale='log')
g.set(xlabel='# of Marxan iterations', ylabel='correlation coefficient')
g.figure.tight_layout()
g.figure.savefig(os.path.join(root,'scripts\dst_04_postprocess', 'corrcoef.svg'))
g.figure.clf()



# add runtime into dataframe and label points with them

# I'm just going to hardcode the numbers in instead of pulling them from
# the text file

runtime = [ # minutes
    0.19,
    0.25,
    0.95,
    7.7,
    59.2,
    67.9
]

df['runtime'] = runtime

# actually, I'll just manually add these in powerpoint
# the text will overlap otherwise