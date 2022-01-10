import main_wrapper as ff
import taxonomic_discordance as td
import utils

if __name__ == '__main__':
    ff.run_procedure(
        model_input_dir="/home/hyperscroll/edwards2016/host/fasta/",
        model_output="/home/hyperscroll/edwards2016/models/wrapped/",
        model_filter="phylum",
        model_minn=6,
        model_maxn=6,
        model_reps=1,
        model_epoch=1,
        model_rm="",
        model_saveVec="",
        general_dim=10,  # 10
        general_length=125,
        general_threads=16,
        workflow_wd="/home/hyperscroll/edwards2016/runs/wrapped_test/",
        workflow_mode="--full",
        workflow_n_vir=200,  # 200
        workflow_n_host=200,  # 200
        workflow_n_nucleotide_threshold=5,
        search_k_nearest=10,  #60
        search_final_rank="wrapped_rank.json",
        bayes_best_score=1000,
        bayes_best_dir="/home/hyperscroll/edwards2016/bayes/"
    )