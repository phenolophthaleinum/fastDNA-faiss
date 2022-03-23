import argparse
import glob
import json
from pathlib import Path, PurePath
import pandas as pd
import re
from collections import defaultdict
# from bokeh.models import BasicTicker, ColorBar, LinearColorMapper, PrintfTickFormatter
# from bokeh.plotting import figure, show
# from bokeh.sampledata.unemployment1948 import data

#fh = open('X:/edwards2016/runs/run-dim_60-len_125-n_600_600-epoch_2-k_7-none_Slim/rank/rank-k_60-new_score.json')
#fh = open('X:/edwards2016/runs/debugF_test/rank/testdebugF_harmonic.json')
#fh = open('bayes_test_avg.json')
p = Path('X:/edwards2016/runs/')
subdirectories = [str(x) for x in p.iterdir() if x.is_dir()]
print(subdirectories)
dirs_to_find = '\n'.join(subdirectories)
dirs = [re.match(r'.*\\debugF_test.*$', x).group(0) for x in subdirectories if re.match(r'.*\\debugF_test.*$', x)]
#filename = 'X:/edwards2016/runs/hybrid_test/rank/testhybrid_avg.json'
print(dirs)

fh = open('virus.json')
vjson = json.load(fh)
fh.close()

fh = open('host.json')
hjson = json.load(fh)
fh.close()

fh = open("hostname.json")
hostname = json.load(fh)
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


final_dict = defaultdict(list)
for d in dirs:
    for filename in glob.glob(f"{d}\\rank\\*.json"):
        key_name = filename.split("\\")[-1]
        print(key_name)
        fh = open(filename)
        db = json.load(fh)
        fh.close()
        d = {}
        d2 = {}
        all_vid = len(list(db.keys()))
        #print(list(db.items())[0])
        for vid, hids in db.items():
            d[vid] = []
            d2[vid] = []
            for taxlevel in TAX_RANKS:
                match = False
                for tup in hids[:1]:
                    hid = tup[0]
                    #mm
                    hid = hostname[hid]['ncbi_id']
                    #print(hid)
                    match = compare_taxranks(vjson[vid], hjson[hid], taxlevel)
                    if match:
                        break
                d[vid].append(match)
                tax_index = hjson[hid]['lineage_ranks'].index(taxlevel)
                d2[vid].append((hjson[hid]['lineage_names'][tax_index], vjson[vid]["host"]["lineage_names"][tax_index], match))
        #print(d2)
        #print(d)
        l = [0 for _i in TAX_RANKS]
        #print(l)
        for vid in d:
            for i, taxlevel in enumerate(TAX_RANKS):
                if d[vid][i]:
                    l[i] += 1
        print(l)


        d = {}
        for i, tax_rank in enumerate(TAX_RANKS):
            #mm
            #d[tax_rank] = l[i] / len(vjson) * 100
            d[tax_rank] = l[i] / all_vid * 100


        print(d)
        final_dict[key_name].append(d)
with open('new_evaluation_final.json', 'w') as fh:
    json.dump(final_dict, fh, indent=4)
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
