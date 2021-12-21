from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
import main_wrapper as ff

# x, y - parameters set to be optimised
pbounds = {'x': (2, 4), 'y': (-3, 3)}

optimizer = BayesianOptimization(
    f=ff.run_procedure(
        model_input_dir="/home/hyperscroll/edwards2016/host/fasta/",
        model_output="/home/hyperscroll/edwards2016/models/wrapped/",
        model_filter="phylum",
        model_minn=6,
        model_maxn=6,
        model_reps=1,
        model_epoch=1,
        model_rm="--rm",
        model_saveVec="",
        general_dim=10,
        general_length=125,
        general_threads=16,
        workflow_wd="/home/hyperscroll/edwards2016/runs/wrapped_test/",
        workflow_mode="--full",
        workflow_n_vir=200,
        workflow_n_host=200,
        workflow_n_nucleotide_threshold=5,
        search_k_nearest=60,
        search_final_rank="wrapped_rank.json"
    ),
    pbounds=pbounds,
    random_state=1,
)
logger = JSONLogger(path="./logs/bayes_opt_log.json")
optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
optimizer.maximize(
    init_points=2,
    n_iter=3,
)
