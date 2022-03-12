import argparse
import json
from pathlib import Path
import pandas as pd
# from bokeh.models import BasicTicker, ColorBar, LinearColorMapper, PrintfTickFormatter
# from bokeh.plotting import figure, show
# from bokeh.sampledata.unemployment1948 import data

# fh = open('X:/edwards2016/runs/run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim/rank/rank-k_60-new_score.json')
fh = open('X:/edwards2016/runs/debugF_test/rank/testdebugF_harmonic.json')
db = json.load(fh)
fh.close()

fh = open('virus.json')
vjson = json.load(fh)
fh.close()

fh = open('host.json')
hjson = json.load(fh)
fh.close()


TAX_RANKS = (
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",   
)


def compare_taxranks(vd, hd, rank):
    vi = vd['host']['lineage_ranks'].index(rank) if rank in vd['host']['lineage_ranks'] else None
    hi = hd['lineage_ranks'].index(rank) if rank in hd['lineage_ranks'] else None
    if vi != None and hi != None:
        if vd['host']['lineage_names'][vi] == hd['lineage_names'][hi]:
            return True
    return False




d = {}
d2 = {}
for vid, hids in db.items():
    d[vid] = []
    d2[vid] = []
    for taxlevel in TAX_RANKS:
        match = False
        for tup in hids[:1]:
            hid = tup[0]
            #print(hid)
            match = compare_taxranks(vjson[vid], hjson[hid], taxlevel)
            if match:
                break
        d[vid].append(match)
        tax_index = hjson[hid]['lineage_ranks'].index(taxlevel)
        d2[vid].append((hjson[hid]['lineage_names'][tax_index], vjson[vid]["host"]["lineage_names"][tax_index], match))
print(d2)
l = [0 for _i in TAX_RANKS]
for vid in d:
    for i, taxlevel in enumerate(TAX_RANKS):
        if d[vid][i]:
            l[i] += 1


d = {}
for i, tax_rank in enumerate(TAX_RANKS):
    d[tax_rank] = l[i] / len(vjson) * 100

print(d)
# print(data)
# data['Year'] = data['Year'].astype(str)
# data = data.set_index('Year')
# data.drop('Annual', axis=1, inplace=True)
# data.columns.name = 'Month'
# print(data)
# df = pd.DataFrame(data.stack(), columns=['rate']).reset_index()
# print(df)
# d2pd = pd.DataFrame.from_dict(d2, orient='index', columns=list(TAX_RANKS))
# d2pd.index.name = "Virus"
# d2pd.columns.name = "Tax level"
# taxranks = list(d2pd.columns)
# viruses = list(d2pd.index)
# print(d2pd)
# TOOLS = "hover,save,pan,box_zoom,reset,wheel_zoom"
# p = figure(title="fastDNA-faiss predictions matches",
#            x_range=viruses, y_range=taxranks,
#            x_axis_location="above", sizing_mode='stretch_both',
#            tools=TOOLS, toolbar_location='below',
#            tooltips=[('predicted', '@Tax level')])
# colors = ["#75968f", "#550b1d"]
# mapper = LinearColorMapper(palette=colors, low=True, high=False)
# p.rect(x="Virus", y="Tax level", width=1, height=1,
#        source=d2pd,
#        fill_color={'field': 'rate', 'transform': mapper},
#        line_color=None)
# show(p)
