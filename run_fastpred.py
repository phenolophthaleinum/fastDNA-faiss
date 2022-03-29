import main_wrapper2 as ff2
import taxonomic_discordance as td
import utils

if __name__ == '__main__':
    ff2.run_procedure(
        model_enable=False,
        model_input_dir="/home/hyperscroll/edwards2016/host/fasta/",
        model_output="/home/hyperscroll/edwards2016/models/debugF_test_Bacillales_bechrun/", # "/home/hyperscroll/edwards2016/models/auto_pred_all/"
        model_filter='order', # species
        model_minn=6,
        model_maxn=6,
        model_reps=1,
        model_epoch=1,
        model_rm="--rm",
        model_saveVec="",
        general_dim=10,  # 10
        general_length=125,
        general_threads=16,
        workflow_wd="/home/hyperscroll/edwards2016/runs/debugF_test_Bacillales_bechrun/", # "/home/hyperscroll/edwards2016/runs/auto_pred_all/"
        workflow_n_vir=200,  # 200
        workflow_n_nucleotide_threshold=5,
        workflow_k_best=10,
        search_final_rank="testdebugF_Bacillales_bechrun.json",
        # search_scoring_func="avg"
    )
