import utils
import os
import glob
import json
from timeit import default_timer as timer
from pathlib import Path

import taxonomic_discordance2 as td
import predict_prob2result as pred

from colorama import Fore, init


# jakie nazwy plikow
# usuwanie niepotrzebnych plikow w dalszej czesci (np. wyrzucenie vectorow oraz sampli hosta, bo i tak jest index)
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
        workflow_n_vir: int,
        workflow_n_nucleotide_threshold: float,
        workflow_k_best: int,
        search_final_rank: str,
        # search_scoring_func: str
        bayes_best_score: float,
        bayes_best_dir: str
):
    # colorama
    init()

    # global timer
    total_start = timer()

    # read config file
    cfg = utils.get_config_obj()

    # clear out models dir
    if len(glob.glob(f"{model_output}*.bin")) != 0:
        os.system(f"rm {model_output}*")
    # run make-model
    # minn = maxn
    print(
        f"python make-model.py -i {model_input_dir} -o {model_output} -f {model_filter} -r {model_reps} -d {general_dim} -l {general_length} --minn {model_maxn} --maxn {model_maxn} -e {model_epoch} -t {general_threads} {model_rm} {model_saveVec}")
    try:
        os.system(
            f"python make-model.py -i {model_input_dir} -o {model_output} -f {model_filter} -r {model_reps} -d {general_dim} -l {general_length} --minn {model_maxn} --maxn {model_maxn} -e {model_epoch} -t {general_threads} {model_rm} {model_saveVec}")
    except RuntimeError:
        exit()
    # assign right model for the next steps
    # also, important, remember to change host and virus values in config.cfg accordingly
    cfg['GENERAL']['active_model'] = glob.glob(f"{model_output}*.bin")[
        0]  # trzeba sie umowic na konwencje nazwywania i organizacji folderow
    with open("config.cfg", "w") as cf:
        cfg.write(cf)
    # run file preprocessing
    try:
        os.system(
            f"python workflow2.py -w {workflow_wd} --length {general_length} --n_vir {workflow_n_vir} -t {general_threads} --n_nucleotide_threshold {workflow_n_nucleotide_threshold} -k {workflow_k_best}")
    except RuntimeError:
        exit()
    # run faiss-search
    # try:
    #     os.system(
    #         f"python predict_prob2result.py -i {workflow_wd}virus/preds/ -o {workflow_wd}rank/{search_final_rank} -s {search_scoring_func}")
    # except RuntimeError:
    #     exit()
    results_table = []

    host_json = Path('host.json')
    with host_json.open() as hj:
        host_dict = json.load(hj)

    hostname_json = Path('hostname.json')
    with hostname_json.open() as hj:
        hostname = json.load(hj)

    virus_json = Path('virus.json')
    with virus_json.open() as hj:
        virus_dict = json.load(hj)
    #
    # # tax_dists = td.load_matrix_parallel('tax_matrix_p2.lzma')
    dists = td.DistanceMatrix(host_dict)

    for fun in pred.scoring:
        preds = pred.score_and_rank(input_dir=f"{workflow_wd}virus/preds/",
                                    scoring_function=fun)
        # accordance = td.taxonomic_accordance_sp(dists, preds, virus_dict)
        score = pred.evaluate_ranks(preds, virus_dict, host_dict, hostname, 'order')
        results_table.append((fun.__name__, preds, score))

    function_ranking = sorted(results_table, key=lambda result: result[2], reverse=True)
    best_func, best_pred, best_score = function_ranking[0]

    total_end = timer()
    total_runtime = total_end - total_start
    print(f"{Fore.GREEN} Total elapsed time: {total_runtime:.6f} seconds {Fore.WHITE}")

    for n, (fun_name, _, score) in enumerate(function_ranking):
        line = f'{n + 1}. {fun_name}: {score}'
        if n:
            print(line)
        else:
            print("\033[1m" + line + " ???" + "\033[0m")
    print(f"[OPT]   Current total best: {bayes_best_score}")
    if best_score > bayes_best_score:
        model_p = Path(f'{bayes_best_dir}{model_output.split("/")[-2]}/')
        if model_p.exists():
            os.system(f'rm -r {str(model_p)}')
        os.system(f"cp -R {workflow_wd} {bayes_best_dir}")
        os.system(f"cp -R {model_output} {bayes_best_dir}")
        dump_path = Path(f'{workflow_wd}rank/{search_final_rank}_{best_func}.json')
        pred.dump_result(output=dump_path,
                         result_rank=best_pred)
        print(f'[OPT]   Top-scoring prediction saved at: {dump_path}')
        # ??
        # pred.evaluate_ranks(best_pred,
        #                     virus_dict,
        #                     host_dict,
        #                     hostname)
    return best_score, best_func
