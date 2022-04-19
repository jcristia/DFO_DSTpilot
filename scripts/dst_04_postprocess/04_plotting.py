# 

import arcpy
import numpy as np
import pandas as pd
import seaborn as sns


root = r'C:\Users\jcristia\Documents\GIS\DFO\DST_pilot'




######################
# targets not met

data =[
    ['prioritzr', 'eco', 0],
    ['prioritzr', 'hu' , 7],
    ['marxan', 'eco', 2],
    ['marxan', 'hu' , 17]
]

df = pd.DataFrame(data, columns = ['pr_ma', 'feature type', 'targets_missed'])

sns.set()
sns.set_style('white')
sns.set_context('notebook', font_scale=1.25)
g = sns.barplot(data=df, x='pr_ma', y='targets_missed', hue='feature type')
g.set(xlabel='Decision support tool', ylabel='# of targets missed')
g.figure.tight_layout()
g.figure.savefig(os.path.join(root,'scripts\dst_04_postprocess', 'targets.svg'))
g.figure.clf()



######################
# % solution in each zone
# (this probably would have been easier to just pull these values manually,
#  oh well)

# hardcode total pus
pus_total = 29231

pr_summ = os.path.join(root, r'scripts\dst_02_prioritizr\outputs\s102_minset_adjtarg\eval_summary.csv')
ma_summ = os.path.join(root, r'scripts\dst_03_marxan\scenario_105_adjtargs_1B\output\output_sum.csv')

best_run = 4 # hardcode the best selected marxan run

pr_summ = pd.read_csv(pr_summ)
ma_summ = pd.read_csv(ma_summ)

df = pd.DataFrame(columns=['DST', 'zone', 'pu_count'])

pr_summ = pr_summ[['zone', 'total_selected']].rename(columns={'total_selected':'pu_count'})
pr_summ = pr_summ[pr_summ.zone != 'overall']
pr_summ['DST'] = 'Prioritizr'
pr_summ.loc[pr_summ.zone == 'zone1', 'zone'] = 'Protection'
pr_summ.loc[pr_summ.zone == 'zone2', 'zone'] = 'Industrial nearshore'
pr_summ.loc[pr_summ.zone == 'zone3', 'zone'] = 'Industrial offshore'
pr_summ.loc[pr_summ.zone == 'zone4', 'zone'] = 'Fishing'
df = df.append(pr_summ)

ma_summ = ma_summ[ma_summ['Run Number']==best_run]
ma_summ = ma_summ.iloc[:, 5:9]
ma_summ = ma_summ.T
ma_summ['zone'] = ma_summ.index
ma_summ = ma_summ.reset_index(drop=True)
ma_summ = ma_summ.rename(columns={3:'pu_count'})
ma_summ['DST'] = 'Marxan'
ma_summ.loc[ma_summ.zone == 'protection PuCount', 'zone'] = 'Protection'
ma_summ.loc[ma_summ.zone == 'industrial_nearshore PuCount', 'zone'] = 'Industrial nearshore'
ma_summ.loc[ma_summ.zone == 'industrial_offshore PuCount', 'zone'] = 'Industrial offshore'
ma_summ.loc[ma_summ.zone == 'fishing PuCount', 'zone'] = 'Fishing'
df = df.append(ma_summ)

df['sol_zone_percent'] = df.pu_count / pus_total * 100

# https://medium.com/dunder-data/automatically-wrap-graph-labels-in-matplotlib-and-seaborn-a48740bc9ce
import textwrap
def wrap_labels(ax, width, break_long_words=False):
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(textwrap.fill(text, width=width,
                      break_long_words=break_long_words))
    ax.set_xticklabels(labels, rotation=0)

sns.set()
sns.set_style('white')
sns.set_context('notebook', font_scale=1.25)
g = sns.barplot(data=df, x='zone', y='sol_zone_percent', hue='DST')
g.set(xlabel='Zone', ylabel='% of total planning units')
wrap_labels(g, 10)
g.figure.savefig(os.path.join(root,'scripts\dst_04_postprocess', 'sol_pus_zone.svg'), bbox_inches="tight", pad_inches=0.25)
g.figure.clf()


# the biggest discrepancy lies in protecting offshore features, which have very 
# linear shapes and might be hard to prioritize correctly

# prioritizr selected more to meet targets. Marxan was unable to do this for
# some reason.
# We could say this is where it needs further calibration - mention feature
# penalty factors.
# So it is likely possible to get them similar, but we wanted to test similarity
# with same input parameters.



######################
# consistency of PUs assigned to each zone
# this will allow us to see if there is a zone type that is assinged differently
# more than others

pr_sol = os.path.join(root, r'scripts\dst_02_prioritizr\outputs\s102_minset_adjtarg\solution.csv')
ma_sol = os.path.join(root, r'scripts\dst_03_marxan\scenario_105_adjtargs_1B\output\output_best.csv')
pr_sol = pd.read_csv(pr_sol)
ma_sol = pd.read_csv(ma_sol)

# reconfigure prioritizr solution
cols = ['uID']
sol_cols = [col for col in pr_sol.columns if col.startswith('solution')]
cols.extend(sol_cols)
pr_sol = pr_sol[cols]
# this is the vectorized numpy way, but the pandas apply/lambda way is
# a bit easier to understand.
conditions = []
outputs = []
i=1
for sol_col in sol_cols:
    condition = pr_sol[sol_col]==1
    conditions.append(condition)
    outputs.append(i)
    i+=1
pr_sol['solution'] = pd.Series(np.select(conditions, outputs, 0))
pr_sol = pr_sol[['uID', 'solution']]

ma_sol['zone'] = ma_sol.zone -1 # adjust marxan zone number to match prioritizr

sol = pr_sol.merge(ma_sol, left_on='uID', right_on='planning_unit')
sol = sol.rename(columns={'solution':'pr_sol', 'zone':'ma_sol', 'uID':'pu'}).drop(columns=['planning_unit'])
sol = sol[sol.pr_sol > 0]

df_summary = pd.DataFrame(columns=['zone', 'total', 'matching'])
for zone in [1,2,3,4]:
    sol_zone = sol[sol.pr_sol==zone]
    sol_zone_total = sol_zone.count()[0]
    sol_zone_match = sol_zone[sol_zone.pr_sol == sol_zone.ma_sol]
    sol_zone_match_total = sol_zone_match.count()[0]
    df2 = {'zone':zone, 'total':sol_zone_total, 'matching':sol_zone_match_total}
    df_summary = df_summary.append(df2, ignore_index=True)

df_summary['perc_matching'] = df_summary.matching / df_summary.total *100

df_summary.loc[df_summary.zone == 1, 'zone'] = 'Protection'
df_summary.loc[df_summary.zone == 2, 'zone'] = 'Industrial nearshore'
df_summary.loc[df_summary.zone == 3, 'zone'] = 'Industrial offshore'
df_summary.loc[df_summary.zone == 4, 'zone'] = 'Fishing'

# https://medium.com/dunder-data/automatically-wrap-graph-labels-in-matplotlib-and-seaborn-a48740bc9ce
import textwrap
def wrap_labels(ax, width, break_long_words=False):
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(textwrap.fill(text, width=width,
                      break_long_words=break_long_words))
    ax.set_xticklabels(labels, rotation=0)

sns.set()
sns.set_style('white')
sns.set_context('notebook', font_scale=1.25)
g = sns.barplot(data=df_summary, x='zone', y='perc_matching')
g.set(xlabel='Zone', ylabel='% PUs matching')
wrap_labels(g, 10)
g.set(ylim=(0, 100))
g.figure.savefig(os.path.join(root,'scripts\dst_04_postprocess', 'pus_matching.svg'), bbox_inches="tight", pad_inches=0.25)
g.figure.clf()




######################
# Edge to area ratio by zone

pr = os.path.join(root, r'spatial\04_outputs\prioritizr.gdb\s102_minset_adjtarg')
ma = os.path.join(root, r'spatial\04_outputs\marxan.gdb\scenario_105_adjtargs_1B')

arcpy.Dissolve_management(pr, 'memory/prdiss', 'solution_1_summary', multi_part='MULTI_PART')
arcpy.Dissolve_management(ma, 'memory/madiss', 'solution_zone', multi_part='MULTI_PART')
arcpy.CalculateGeometryAttributes_management('memory/prdiss', [['area', 'AREA'], ['length', 'PERIMETER_LENGTH']])
arcpy.CalculateGeometryAttributes_management('memory/madiss', [['area', 'AREA'], ['length', 'PERIMETER_LENGTH']])

field_names = [i.name for i in arcpy.ListFields('memory/prdiss') if i.type not in ['OID', 'Geometry']]
cursor = arcpy.da.SearchCursor('memory/prdiss', field_names)
df_pr = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df_pr = df_pr.rename(columns={'solution_1_summary':'zone'})
df_pr.loc[df_pr.zone == 1, 'zone'] = 'Protection'
df_pr.loc[df_pr.zone == 2, 'zone'] = 'Industrial nearshore'
df_pr.loc[df_pr.zone == 3, 'zone'] = 'Industrial offshore'
df_pr.loc[df_pr.zone == 4, 'zone'] = 'Fishing'
df_pr = df_pr[df_pr.zone != 0]
df_pr['length'] = df_pr.length / 1000
df_pr['area'] = df_pr.area / 1000000
df_total = {'zone':'Overall', 'area':df_pr.area.sum(), 'length':df_pr.length.sum()}
df_pr = df_pr.append(df_total, ignore_index=True)
df_pr['DST'] = 'Prioritizr'


field_names = [i.name for i in arcpy.ListFields('memory/madiss') if i.type not in ['OID', 'Geometry']]
cursor = arcpy.da.SearchCursor('memory/madiss', field_names)
df_ma = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df_ma = df_ma.rename(columns={'solution_zone':'zone'})
df_ma.loc[df_ma.zone == 1, 'zone'] = 'Protection'
df_ma.loc[df_ma.zone == 2, 'zone'] = 'Industrial nearshore'
df_ma.loc[df_ma.zone == 3, 'zone'] = 'Industrial offshore'
df_ma.loc[df_ma.zone == 4, 'zone'] = 'Fishing'
df_ma = df_ma[df_ma.zone != 0]
df_ma['length'] = df_ma.length / 1000
df_ma['area'] = df_ma.area / 1000000
df_total = {'zone':'Overall', 'area':df_ma.area.sum(), 'length':df_ma.length.sum()}
df_ma = df_ma.append(df_total, ignore_index=True)
df_ma['DST'] = 'Marxan'

df_all = df_pr.append(df_ma)
df_all['ratio'] = df_all.length / df_all.area

import textwrap
def wrap_labels(ax, width, break_long_words=False):
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(textwrap.fill(text, width=width,
                      break_long_words=break_long_words))
    ax.set_xticklabels(labels, rotation=0)

sns.set()
sns.set_style('white')
sns.set_context('notebook', font_scale=1.25)
g = sns.barplot(data=df_all, x='zone', y='ratio', hue='DST')
g.set(xlabel='Zone', ylabel='Edge to area ratio')
wrap_labels(g, 10)
g.figure.savefig(os.path.join(root,'scripts\dst_04_postprocess', 'edge_area.svg'), bbox_inches="tight", pad_inches=0.25)
g.figure.clf()


# see if any of the length match (or are a factor of) the boundary numbers in
# the result csvs. It's not clear to me.
#
# In the prioritizr eval_boundary.csv, the total length is very close, but not
# exactly the same. I wonder if this is with how donut holes are summed.
# In the Marxan output_sum.csv Connectivity Strength, it doesn't seem to match
# at all unless it somehow knew to convert to kilometers.



# Add the prioritizr with boundary penalty scenario to the dataframe
prb = os.path.join(root, r'spatial\04_outputs\prioritizr.gdb\s103_minset_adjtarg_bp')
arcpy.Dissolve_management(prb, 'memory/prbdiss', 'solution_1_summary', multi_part='MULTI_PART')
arcpy.CalculateGeometryAttributes_management('memory/prbdiss', [['area', 'AREA'], ['length', 'PERIMETER_LENGTH']])

field_names = [i.name for i in arcpy.ListFields('memory/prbdiss') if i.type not in ['OID', 'Geometry']]
cursor = arcpy.da.SearchCursor('memory/prbdiss', field_names)
df_prb = pd.DataFrame(data=[row for row in cursor], columns=field_names)
df_prb = df_prb.rename(columns={'solution_1_summary':'zone'})
df_prb.loc[df_prb.zone == 1, 'zone'] = 'Protection'
df_prb.loc[df_prb.zone == 2, 'zone'] = 'Industrial nearshore'
df_prb.loc[df_prb.zone == 3, 'zone'] = 'Industrial offshore'
df_prb.loc[df_prb.zone == 4, 'zone'] = 'Fishing'
df_prb = df_prb[df_prb.zone != 0]
df_prb['length'] = df_prb.length / 1000
df_prb['area'] = df_prb.area / 1000000
df_total = {'zone':'Overall', 'area':df_prb.area.sum(), 'length':df_prb.length.sum()}
df_prb = df_prb.append(df_total, ignore_index=True)
df_prb['DST'] = 'Prioritizr BP'
df_all2 = df_all.append(df_prb)
df_all2['ratio'] = df_all2.length / df_all2.area

sns.set()
sns.set_style('white')
sns.set_context('notebook', font_scale=1.25)
g = sns.barplot(data=df_all2, x='zone', y='ratio', hue='DST')
g.set(xlabel='Zone', ylabel='Edge to area ratio')
wrap_labels(g, 10)
g.figure.savefig(os.path.join(root,'scripts\dst_04_postprocess', 'edge_area_BP.svg'), bbox_inches="tight", pad_inches=0.25)
g.figure.clf()
