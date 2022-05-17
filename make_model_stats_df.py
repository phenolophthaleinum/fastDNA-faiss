import pandas as pd
import utils
from Bio import SeqIO
from timeit import default_timer as timer

start = timer()
hostname = utils.get_hostname_data()
filename = "X:/edwards2016/models/stats_test/random_labels-species-dim_10-len_125.txt"


def get_genome_length(ncbi_id):
    record, = SeqIO.parse(f"X:/edwards2016/host/fasta/{ncbi_id}.fna", 'fasta')
    return len(record.seq)


with open(filename) as file:
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]
df = pd.DataFrame()
df['organism'] = lines
df['phylum'] = df['organism'].apply(lambda x: hostname[x]['lineage_names'][1])
df['class'] = df['organism'].apply(lambda x: hostname[x]['lineage_names'][2])
df['order'] = df['organism'].apply(lambda x: hostname[x]['lineage_names'][3])
df['family'] = df['organism'].apply(lambda x: hostname[x]['lineage_names'][4])
df['genus'] = df['organism'].apply(lambda x: hostname[x]['lineage_names'][5])
df['genome length'] = df['organism'].apply(lambda x: get_genome_length(hostname[x]['ncbi_id']))

df_genus = pd.DataFrame()
sub_df = df[df['genus'] == '']
df_genus_unnamed = sub_df['organism']
df_genus['genome counts'] = df['genus'].value_counts()
df_genus['genome length avg'] = df.groupby('genus').mean()
df_genus['genome length std'] = df.groupby('genus').std()
df_genus['genome length std'].fillna(0, inplace=True)
df_genus['taxlevel'] = 'genus'
df_genus.index.name = 'names'
df_genus.reset_index(inplace=True)
# df_genus = df_genus.reset_index()

df_family = pd.DataFrame()
sub_df = df[df['family'] == '']
df_family_unnamed = sub_df['organism']
df_family['genome counts'] = df['family'].value_counts()
df_family['genome length avg'] = df.groupby('family').mean()
df_family['genome length std'] = df.groupby('family').std()
df_family['genome length std'].fillna(0, inplace=True)
df_family['taxlevel'] = 'family'
df_family.index.name = 'names'
df_family.reset_index(inplace=True)
# df_family = df_family.reset_index()

df_order = pd.DataFrame()
sub_df = df[df['order'] == '']
df_order_unnamed = sub_df['organism']
df_order['genome counts'] = df['order'].value_counts()
df_order['genome length avg'] = df.groupby('order').mean()
df_order['genome length std'] = df.groupby('order').std()
df_order['genome length std'].fillna(0, inplace=True)
df_order['taxlevel'] = 'order'
df_order.index.name = 'names'
df_order.reset_index(inplace=True)
# df_order = df_order.reset_index()

df_class = pd.DataFrame()
sub_df = df[df['class'] == '']
df_class_unnamed = sub_df['organism']
df_class['genome counts'] = df['class'].value_counts()
df_class['genome length avg'] = df.groupby('class').mean()
df_class['genome length std'] = df.groupby('class').std()
df_class['genome length std'].fillna(0, inplace=True)
df_class['taxlevel'] = 'class'
df_class.index.name = 'names'
df_class.reset_index(inplace=True)
# df_class = df_class.reset_index()

df_phylum = pd.DataFrame()
sub_df = df[df['phylum'] == '']
df_phylum_unnamed = sub_df['organism']
df_phylum['genome counts'] = df['phylum'].value_counts()
df_phylum['genome length avg'] = df.groupby('phylum').mean()
df_phylum['genome length std'] = df.groupby('phylum').std()
df_phylum['genome length std'].fillna(0, inplace=True)
df_phylum['taxlevel'] = 'phylum'
df_phylum.index.name = 'names'
df_phylum.reset_index(inplace=True)
# df_phylum = df_phylum.reset_index()

df_final = pd.concat([df_genus, df_family, df_order, df_class, df_phylum])

end = timer()
runtime = end - start
print(f"Done in {runtime:.6f}")
