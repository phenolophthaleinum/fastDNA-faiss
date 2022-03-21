import os
import utils

virus_data = utils.get_virus_data()
config = utils.get_config()
filter = 'Clostridiales'
target = f'/home/hyperscroll/edwards2016/virus/subsets/{filter.lower()}/'
os.system(f'mkdir {target}')
print(target)
for virus in virus_data:
    if filter in virus_data[virus]["host"]["lineage_names"]:
        os.system(f'cp {config["VIRUS"]["virus_genomes"]}{virus}.fna {target}')
