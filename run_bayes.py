from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
import main_wrapper as ff

best = 0


def bayes():
    global best

    def procedure_wrapper(model_maxn, general_dim, general_length):
        model_maxn_d = int(model_maxn)
        general_dim_d = int(general_dim)
        general_length_d = int(general_length)
        global best
        ff_score = ff.run_procedure(
            model_input_dir="/home/hyperscroll/edwards2016/host/fasta/",  # slim
            model_output="/home/hyperscroll/edwards2016/models/wrapped/",
            model_filter="phylum",  # none
            model_minn=6,
            model_maxn=model_maxn_d,
            model_reps=1,
            model_epoch=1,  # 10
            model_rm="--rm",
            model_saveVec="",
            general_dim=general_dim_d,
            general_length=general_length_d,
            general_threads=16,
            workflow_wd="/home/hyperscroll/edwards2016/runs/wrapped_test/",
            workflow_mode="--full",
            workflow_n_vir=100,  # 500
            workflow_n_host=100,  # 500
            workflow_n_nucleotide_threshold=5,
            search_k_nearest=10,  # 100
            search_final_rank="wrapped_rank.json",
            bayes_best_score=best,
            bayes_best_dir="/home/hyperscroll/edwards2016/bayes/"
        )
        if ff_score > best:
            best = ff_score
            print(f"[OPT]   New best: {best}")
        return ff_score

    # x, y - parameters set to be optimised
    pbounds = {'model_maxn': (6, 7),  # (6, 10)
               'general_dim': (10, 12),  # (20, 200)
               'general_length': (50, 100)}  # (50, 250)

    optimizer = BayesianOptimization(
        f=procedure_wrapper,
        pbounds=pbounds,
        random_state=1,
        verbose=2
    )
    logger = JSONLogger(path="./logs/bayes_opt_log.json")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
    optimizer.maximize(
        init_points=2,  # 30
        n_iter=3,  # 90
    )

    return optimizer.max


if __name__ == "__main__":
    print(bayes())
