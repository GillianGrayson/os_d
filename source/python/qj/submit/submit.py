import pathlib
import os.path
from config.config import *
from config.run import *
from config.details import *
from config.random import *
from config.params import *
from experiments.times_between_jumps import *

run = Run(quantum_system=QuantumSystem.photonic,
          process=Process.evolution_std,
          propagation=Propagation.quantum_jumps)

details = Details(step_metrics=16,
                  num_periods_trans=1000,
                  num_periods_obser=10000)
if run.propagation is Propagation.quantum_jumps:
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
params = Params(params_size,params_dissipation, params_driving, params_init_cond, params_model)

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


                diss_gamma_str = str(format(diss_gamma, '0.4f'))
                diss_phase_str = str(format(diss_phase, '0.4f'))

                dimer_drv_ampl_str = str(format(dimer_drv_ampl, '0.4f'))
                dimer_drv_freq_str = str(format(dimer_drv_freq, '0.4f'))
                dimer_drv_phase_str = str(format(dimer_drv_phase, '0.4f'))
                dimer_prm_E_str = str(format(dimer_prm_E, '0.4f'))
                dimer_prm_U_str = str(format(dimer_prm_U, '0.4f'))
                dimer_prm_J_str = str(format(dimer_prm_J, '0.4f'))

                jcs_drv_part_1_str = str(format(jcs_drv_part_1, '0.4f'))
                jcs_drv_part_2_str = str(format(jcs_drv_part_2, '0.4f'))
                jcs_drv_ampl_str = str(format(jcs_drv_ampl, '0.4f'))
                jcs_prm_alpha_str = str(format(jcs_prm_alpha, '0.4f'))

                start_seed = 0
                finish_seed = num_runs * num_trajectories
                step_seed = num_trajectories

                local_path = ''
                file_suffix = ''
                if sys_id == 0:

                    local_path = \
                        '/main_' + str(sys_id) + '_' + str(task_id) + '_' + str(prop_id) + \
                        '/run_' + str(ex_deep) + '_' + str(rk_ns) + '_' + str(num_tp_periods) + '_' + str(
                            num_obs_periods) + \
                        '/N_' + str(N) + \
                        '/diss_' + str(diss_type) + '_' + diss_gamma_str + '_' + diss_phase_str + \
                        '/drv_' + str(
                            dimer_drv_type) + '_' + dimer_drv_ampl_str + '_' + dimer_drv_freq_str + '_' + dimer_drv_phase_str + \
                        '/prm_' + dimer_prm_E_str + '_' + dimer_prm_U_str + '_' + dimer_prm_J_str + \
                        '/start_' + str(start_type) + '_' + str(start_state)

                elif sys_id == 1:

                    local_path = \
                        '/main_' + str(sys_id) + '_' + str(task_id) + '_' + str(prop_id) + \
                        '/run_' + str(ex_deep) + '_' + str(rk_ns) + '_' + str(num_tp_periods) + '_' + str(
                            num_obs_periods) + \
                        '/N_' + str(N) + \
                        '/diss_' + str(diss_type) + '_' + diss_gamma_str + '_' + diss_phase_str + \
                        '/drv_' + jcs_drv_part_1_str + '_' + jcs_drv_part_2_str + '_' + jcs_drv_ampl_str + \
                        '/prm_' + jcs_prm_alpha_str + \
                        '/start_' + str(start_type) + '_' + str(start_state)

                for ss in range(start_seed, step_seed, finish_seed):
                    fn_path = get_dir(local_path, ss, type)
                    pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

                    file_config = open(fn_path + '/config.txt', 'w')

                    file_config.write('sys_id ' + str(sys_id) + '\n')
                    file_config.write('task_id ' + str(task_id) + '\n')
                    file_config.write('prop_id ' + str(prop_id) + '\n')
                    file_config.write('is_debug ' + str(is_debug) + '\n')
                    file_config.write('is_pp ' + str(is_pp) + '\n')
                    file_config.write('init_fn ' + init_fn + '\n')
                    file_config.write('path ' + path + '\n')
                    file_config.write('seed ' + str(ss) + '\n')
                    file_config.write('mns ' + str(mns) + '\n')
                    file_config.write('num_threads ' + str(num_threads) + '\n')
                    file_config.write('num_trajectories ' + str(num_trajectories) + '\n')
                    file_config.write('num_tp_periods ' + str(num_tp_periods) + '\n')
                    file_config.write('num_obs_periods ' + str(num_obs_periods) + '\n')
                    file_config.write('ex_deep ' + str(ex_deep) + '\n')
                    file_config.write('rk_ns ' + str(rk_ns))

                    file_config.close()

                    file_params = open(fn_path + '/params.txt', 'w')

                    file_params.write('lpn_type ' + str(lpn_type) + '\n')
                    file_params.write('lpn_eps ' + str(lpn_eps) + '\n')
                    file_params.write('lpn_eps_change ' + str(lpn_eps_change) + '\n')
                    file_params.write('lpn_delta_s_high ' + str(lpn_delta_s_high) + '\n')
                    file_params.write('lpn_delta_s_low ' + str(lpn_delta_s_low) + '\n')
                    file_params.write('lpn_delta_f_high ' + str(lpn_delta_f_high) + '\n')
                    file_params.write('lpn_delta_f_low ' + str(lpn_delta_f_low) + '\n')
                    file_params.write('dump_obs ' + str(dump_obs) + '\n')
                    file_params.write('dump_adr_sep ' + str(dump_adr_sep) + '\n')
                    file_params.write('dump_adr_avg ' + str(dump_adr_avg) + '\n')
                    file_params.write('dump_evo_sep ' + str(dump_evo_sep) + '\n')
                    file_params.write('dump_evo_avg ' + str(dump_evo_avg) + '\n')
                    file_params.write('dump_type ' + str(dump_type) + '\n')
                    file_params.write('dump_num ' + str(dump_num) + '\n')
                    file_params.write('N ' + str(N) + '\n')
                    file_params.write('diss_type ' + str(diss_type) + '\n')
                    file_params.write('diss_gamma ' + str(diss_gamma) + '\n')
                    file_params.write('diss_phase ' + str(diss_phase) + '\n')
                    file_params.write('dimer_drv_type ' + str(dimer_drv_type) + '\n')
                    file_params.write('dimer_drv_ampl ' + str(dimer_drv_ampl) + '\n')
                    file_params.write('dimer_drv_freq ' + str(dimer_drv_freq) + '\n')
                    file_params.write('dimer_drv_phase ' + str(dimer_drv_phase) + '\n')
                    file_params.write('dimer_prm_E ' + str(dimer_prm_E) + '\n')
                    file_params.write('dimer_prm_U ' + str(dimer_prm_U) + '\n')
                    file_params.write('dimer_prm_J ' + str(dimer_prm_J) + '\n')
                    file_params.write('jcs_drv_part_1 ' + str(jcs_drv_part_1) + '\n')
                    file_params.write('jcs_drv_part_2 ' + str(jcs_drv_part_2) + '\n')
                    file_params.write('jcs_drv_ampl ' + str(jcs_drv_ampl) + '\n')
                    file_params.write('jcs_prm_alpha ' + str(jcs_prm_alpha) + '\n')
                    file_params.write('start_type ' + str(start_type) + '\n')
                    file_params.write('start_state ' + str(start_state) + '\n')
                    file_params.write('cd_dim ' + str(cd_dim) + '\n')
                    file_params.write('cd_eps ' + str(cd_eps) + '\n')
                    file_params.write('deep_num_steps ' + str(deep_num_steps) + '\n')
                    file_params.write('deep_dump ' + str(deep_dump) + '\n')
                    file_params.write('jump ' + str(jump))

                    file_params.close()

                    fn_suffix = ''
                    if sys_id == 0:

                        fn_suffix = \
                            'config(' + str(sys_id) + '_' + str(task_id) + '_' + str(prop_id) + ')_' + \
                            'rnd(' + str(ss) + '_' + str(mns) + ')_' + \
                            'N(' + str(N) + ')_' + \
                            'diss(' + str(diss_type) + '_' + diss_gamma_str + '_' + diss_phase_str + ')_' + \
                            'drv(' + str(
                                dimer_drv_type) + '_' + dimer_drv_ampl_str + '_' + dimer_drv_freq_str + '_' + dimer_drv_phase_str + ')_' + \
                            'prm(' + dimer_prm_E_str + '_' + dimer_prm_U_str + '_' + dimer_prm_J_str + ')_' + \
                            'start(' + str(start_type) + '_' + str(start_state) + ')'

                    elif sys_id == 1:

                        fn_suffix = \
                            'config(' + str(sys_id) + '_' + str(task_id) + '_' + str(prop_id) + ')_' + \
                            'rnd(' + str(ss) + '_' + str(mns) + ')_' + \
                            'N(' + str(N) + ')_' + \
                            'diss(' + str(diss_type) + '_' + diss_gamma_str + '_' + diss_phase_str + ')_' + \
                            'drv(' + jcs_drv_part_1_str + '_' + jcs_drv_part_2_str + '_' + jcs_drv_ampl_str + ')_' + \
                            'prm(' + jcs_prm_alpha_str + ')_' + \
                            'start(' + str(start_type) + '_' + str(start_state) + ')'

                    fn_test = ''
                    if sys_id == 0:
                        fn_test = fn_path + '/mean_' + fn_suffix + '.txt'
                    elif sys_id == 1:
                        fn_test = fn_path + '/spec_' + fn_suffix + '.txt'

                    if not os.path.isfile(fn_test):
                        os.system('sbatch run_unn.sh ' + fn_path)








