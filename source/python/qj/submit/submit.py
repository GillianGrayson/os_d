import pathlib
from config.config import *
from config.run import *
from config.details import *
from config.random import *
from config.params import *
from experiments.times_between_jumps import *
from infrastructure.path import *
from infrastructure.save import *
import os

run = Run(quantum_system=QuantumSystem.photonic,
          process=Process.evolution_std,
          propagation=Propagation.quantum_jumps)

details = Details(step_metrics=16,
                  num_periods_trans=1000,
                  num_periods_obser=10000)
if run.propagation is Propagation.runge_kutta_4:
    details.step_metrics = 10000

random = Random(seed=0,
                num_seeds=1000000)

params_size = {'N': 200}
params_dissipation = {'diss_type': 0,
                      'diss_gamma': 0.1,
                      'diss_phase': 0.0}
params_driving = {'jcs_drv_part_1': 0.98,
                  'jcs_drv_part_2': 1.00,
                  'jcs_drv_ampl': 3.2}
params_init_cond = {'start_type': 0,
                    'start_state': 0}
params_model = {'jcs_prm_alpha': 5.0}
if run.quantum_system is QuantumSystem.photonic:
    params_driving = {'dimer_drv_type': 0,
                      'dimer_drv_ampl': 1.50,
                      'dimer_drv_freq': 1.0,
                      'dimer_drv_phase': 0.0}
    params_model = {'dimer_prm_E': 1.0,
                    'dimer_prm_U': 0.5,
                    'dimer_prm_J': 1.0}
params = Params(params_size, params_dissipation, params_driving, params_model, params_init_cond)

auxiliary = {'is_debug': 0,
             'is_pp': 0,
             'init_fn': '',
             'path': '',
             'num_threads': 32,
             'num_trajectories': 32,
             'lpn_type': 0,
             'lpn_eps_deep': 10,
             'lpn_eps_high': 0,
             'lpn_eps_low': -4,
             'lpn_delta_s_high': 1.0e-3,
             'lpn_delta_s_low': 1.0e-4,
             'lpn_delta_f_high': 1.0e-2,
             'lpn_delta_f_low': 1.0e-10,
             'dump_obs': 1,
             'dump_adr_sep': 0,
             'dump_adr_avg': 0,
             'dump_evo_sep': 1,
             'dump_evo_avg': 0,
             'dump_type': 0,
             'dump_num': 1,
             'cd_dim': 1,
             'cd_eps': 1.0e-4,
             'deep_num_steps': 10000,
             'jump': 2,
             'jumps_counts': 1000}

config = Config(run, details, random, params, auxiliary)

[xs, ys, configs] = get_space(config)

for curr_config in configs:
    path = get_data_path(curr_config)
    suffix = 'rnd(' + str(curr_config.random.seed) + '_' + str(curr_config.random.num_seeds) + ').txt'
    for fn in os.listdir(path):
        if fn.endswith(suffix):

            create_file(curr_config)
            os.system('sbatch run_unn.sh ' + path)

            break