from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
import main_wrapper2_bay as ff
import os
from bayes_opt.observer import _Tracker
import json


ff_func = ""
best = 0


class FastpredLogger(_Tracker):
    def __init__(self, path, reset=True):
        self._path = path if path[-5:] == ".json" else path + ".json"
        if reset:
            try:
                os.remove(self._path)
            except OSError:
                pass
        super(FastpredLogger, self).__init__()

    def update(self, event, instance):
        data = dict(instance.res[-1])

        now, time_elapsed, time_delta = self._time_metrics()
        data["datetime"] = {
            "datetime": now,
            "elapsed": time_elapsed,
            "delta": time_delta,
        }
        data["function"] = ff_func

        # print(dict(data))

        with open(self._path, "a") as f:
            f.write(json.dumps(dict(data)) + "\n")

        self._update_tracker(event, instance)


def bayes():
    global best
    global ff_func

    def procedure_wrapper(model_maxn, general_dim, general_length):
        model_maxn_d = int(model_maxn)
        general_dim_d = int(general_dim)
        general_length_d = int(general_length)
        global best
        global ff_func
        ff_score, ff_func = ff.run_procedure(
            model_input_dir="/home/hyperscroll/edwards2016/host/fasta/",  # slim
            model_output="/home/hyperscroll/edwards2016/models/bayes_test_rework/",
            model_filter="order",  # none
            model_minn=6,
            model_maxn=model_maxn_d,
            model_reps=1,
            model_epoch=1,  # 10
            model_rm="--rm",
            model_saveVec="",
            general_dim=general_dim_d,
            general_length=general_length_d,
            general_threads=16,
            workflow_wd="/home/hyperscroll/edwards2016/runs/bayes_test_rework/",
            workflow_n_vir=100,  # 500
            workflow_n_nucleotide_threshold=5,
            workflow_k_best=10,
            search_final_rank="bayes_test_rework.json",
            bayes_best_score=best,
            bayes_best_dir="/home/hyperscroll/edwards2016/bayes_rework/"
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
    # logger = JSONLogger(path="./logs/bayes_opt_log.json")
    observer = FastpredLogger(path="./logs/bayes_opt_log_rework.json")
    # optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
    optimizer.subscribe(Events.OPTIMIZATION_STEP, observer)
    optimizer.maximize(
        init_points=2,  # 30
        n_iter=3,  # 90
    )

    return optimizer.max


if __name__ == "__main__":
    print(bayes())
