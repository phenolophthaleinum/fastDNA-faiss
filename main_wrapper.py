import argparse
import configparser
import utils
import os
import glob
from timeit import default_timer as timer

from colorama import Fore, init


# jakie nazwy plikow
# usuwanie niepotrzebnych plikow w dalszej czesci (np. wyrzucenie vectorow oraz sampli hosta, bo i tak jest index)
# usuwanie N w sekwencjach
def run_procedure(
        model_input_dir: str,
        model_output: str,
        model_filter: str,
        model_reps: int,
        model_minn: int,
        model_maxn: int,
        model_epoch: int,
        model_rm: str,
        model_saveVec: str,
        general_dim: int,
        general_length: int,
        general_threads: int,
        workflow_wd: str,
        workflow_mode: str,
        workflow_n_vir: int,
        workflow_n_host: int,
        search_k_nearest: int,
        search_final_rank: str
):
    # colorama
    init()

    # global timer
    total_start = timer()

    # read config file
    cfg = utils.get_config_obj()
    # run make-model
    print(f"python make-model.py -i {model_input_dir} -o {model_output} -f {model_filter} -r {model_reps} -d {general_dim} -l {general_length} --minn {model_minn} --maxn {model_maxn} -e {model_epoch} -t {general_threads} {model_rm} {model_saveVec}")
    try:
        os.system(f"python make-model.py -i {model_input_dir} -o {model_output} -f {model_filter} -r {model_reps} -d {general_dim} -l {general_length} --minn {model_minn} --maxn {model_maxn} -e {model_epoch} -t {general_threads} {model_rm} {model_saveVec}")
    except RuntimeError:
        exit()
    # assign right model for the next steps
    # also, important, remember to change host and virus values in config.cfg accordingly
    cfg['GENERAL']['active_model'] = glob.glob(f"{model_output}*.bin")[0] # trzeba sie umowic na konwencje nazwywania i organizacji folderow
    with open("config.cfg", "w") as cf:
        cfg.write(cf)
    # run file preprocessing
    try:
        os.system(f"python workflow.py -w {workflow_wd} {workflow_mode} --length {general_length} -d {general_dim} --n_vir {workflow_n_vir} --n_host {workflow_n_host} -t {general_threads}")
    except RuntimeError:
        exit()
    # run faiss-search
    try:
        os.system(f"python faiss_search.py -i {workflow_wd}virus/vectors/ -o {workflow_wd}rank/{search_final_rank} -d {general_dim} --n_samples {workflow_n_host} -k {search_k_nearest} -f {glob.glob(f'{workflow_wd}host/index/*.index')[0]} -m {glob.glob(f'{workflow_wd}host/maps/*.json')[0]}")
    except RuntimeError:
        exit()

    # run evaluation
    # ?

    total_end = timer()
    total_runtime = total_end - total_start
    print(f"{Fore.GREEN} Total elapsed time: {total_runtime:.6f} seconds")
