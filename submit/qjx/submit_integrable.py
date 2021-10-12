import pathlib
from Infrastructure.file_system import *
import os.path
import numpy as np

type = FSType.mpipks_sd

num_runs = 1

medium = 0

seed_start = 1
seed_shift = 1
seed_num = 50

Ts = [1.0, 2.0, 5.0, 10.0]
deltas = [1e-6, 1e-5, 1e-4, 1e-3]

for T in Ts:
    for delta in deltas:
        for seed_id in range(0, seed_num):

            integrable_seed = seed_start + seed_id * seed_shift

            print(f"T: {T}")
            print(f"delta: {delta}")
            print(f"integrable_seed: {integrable_seed}")

            sys_id = 5
            task_id = 7
            prop_id = 0
            is_debug = 0
            is_pp = 0
            init_fn = ''
            path = ''
            seed = 0
            mns = 1000000
            num_threads = 1
            num_trajectories = 100
            num_tp_periods = 100
            num_obs_periods = 100
            ex_deep = 16
            rk_ns = 10000

            num_random_obs = 1
            random_obs_seed = integrable_seed
            random_obs_mns = 1000000
            random_obs_type = 2
            lpn_type = 0
            lpn_eps_deep = 100
            lpn_eps_error = 1.0e-10
            lpn_eps_high = 10
            lpn_eps_low = -10
            lpn_delta_s = delta
            lpn_delta_f_high = delta
            lpn_delta_f_low = delta
            lambda_per_periods = 1
            save_lambdas = 0
            num_lambdas_periods = 2
            dump_obs = 1
            dump_phi = 0
            dump_phi_evo = 0
            dump_adr_sep = 0
            dump_adr_avg = 0
            dump_evo_sep = 1
            dump_evo_avg = 0
            dump_type = 0
            dump_num = 1
            diss_type = 1
            diss_gamma = 0.1
            diss_phase = 0.0
            integrable_N = 7
            integrable_seed = integrable_seed
            integrable_num_seeds = 1000000
            integrable_tau = 1
            integrable_k = -1
            integrable_T = T
            start_type = 0
            start_state = 0
            cd_dim = 1
            cd_eps = 1.0e-8
            deep_num_steps = 1000
            jump = 0
            jumps_counts = 0

            start_seed = 0
            finish_seed = num_runs * num_trajectories
            step_seed = num_trajectories

            local_path = f"/main_{sys_id}_{task_id}_{prop_id}/T_{integrable_T:0.4f}"

            if task_id == 7:
                local_path += f"/lpn_{lpn_type}_{np.log10(lpn_delta_s):0.4f}"

            local_path += f"/run_{ex_deep}_{rk_ns}_{num_tp_periods}_{num_obs_periods}" + \
                          f"/N_{integrable_N}_tau_{integrable_tau:d}_k_{integrable_k:d}" + \
                          f'/seed_{integrable_seed}'

            for ss in range(start_seed, finish_seed, step_seed):
                print("ss = " + str(ss))
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
                file_params.write('num_random_obs ' + str(num_random_obs) + '\n')
                file_params.write('random_obs_seed ' + str(random_obs_seed) + '\n')
                file_params.write('random_obs_mns ' + str(random_obs_mns) + '\n')
                file_params.write('random_obs_type ' + str(random_obs_type) + '\n')
                file_params.write('lpn_type ' + str(lpn_type) + '\n')
                file_params.write('lpn_eps_deep ' + str(lpn_eps_deep) + '\n')
                file_params.write('lpn_eps_error ' + str(lpn_eps_error) + '\n')
                file_params.write('lpn_eps_high ' + str(lpn_eps_high) + '\n')
                file_params.write('lpn_eps_low ' + str(lpn_eps_low) + '\n')
                file_params.write('lpn_delta_s ' + str(lpn_delta_s) + '\n')
                file_params.write('lpn_delta_f_high ' + str(lpn_delta_f_high) + '\n')
                file_params.write('lpn_delta_f_low ' + str(lpn_delta_f_low) + '\n')
                file_params.write('lambda_per_periods ' + str(lambda_per_periods) + '\n')
                file_params.write('save_lambdas ' + str(save_lambdas) + '\n')
                file_params.write('num_lambdas_periods ' + str(num_lambdas_periods) + '\n')
                file_params.write('dump_obs ' + str(dump_obs) + '\n')
                file_params.write('dump_phi ' + str(dump_phi) + '\n')
                file_params.write('dump_phi_evo ' + str(dump_phi_evo) + '\n')
                file_params.write('dump_adr_sep ' + str(dump_adr_sep) + '\n')
                file_params.write('dump_adr_avg ' + str(dump_adr_avg) + '\n')
                file_params.write('dump_evo_sep ' + str(dump_evo_sep) + '\n')
                file_params.write('dump_evo_avg ' + str(dump_evo_avg) + '\n')
                file_params.write('dump_type ' + str(dump_type) + '\n')
                file_params.write('dump_num ' + str(dump_num) + '\n')
                file_params.write('diss_type ' + str(diss_type) + '\n')
                file_params.write('diss_gamma ' + str(diss_gamma) + '\n')
                file_params.write('diss_phase ' + str(diss_phase) + '\n')
                file_params.write('integrable_N ' + str(integrable_N) + '\n')
                file_params.write('integrable_seed ' + str(integrable_seed) + '\n')
                file_params.write('integrable_num_seeds ' + str(integrable_num_seeds) + '\n')
                file_params.write('integrable_tau ' + str(integrable_tau) + '\n')
                file_params.write('integrable_k ' + str(integrable_k) + '\n')
                file_params.write('integrable_T ' + str(integrable_T) + '\n')
                file_params.write('start_type ' + str(start_type) + '\n')
                file_params.write('start_state ' + str(start_state) + '\n')
                file_params.write('cd_dim ' + str(cd_dim) + '\n')
                file_params.write('cd_eps ' + str(cd_eps) + '\n')
                file_params.write('deep_num_steps ' + str(deep_num_steps) + '\n')
                file_params.write('jump ' + str(jump) + '\n')
                file_params.write('jumps_counts ' + str(jumps_counts))

                file_params.close()

                fn_suffix = f"setup({sys_id}_{task_id}_{prop_id})_rnd({ss}_{mns})_N({integrable_N})_seed({integrable_seed})_tau({integrable_tau:d})_k({integrable_k:d})_T({integrable_T:0.4f})"
                if task_id == 7:
                    fn_suffix += f"_lpn({lpn_type}_{np.log10(lpn_delta_s):0.4f}_{np.log10(lpn_delta_f_high):0.4f}_{np.log10(lpn_delta_f_low):0.4f})"

                fn_test = fn_path + '/spec_' + fn_suffix + '.txt'

                if not os.path.isfile(fn_test):
                    print(fn_test)
                    if type == FSType.cluster:
                        os.system('sbatch run_unn.sh ' + fn_path)
                    elif type == FSType.mpipks_sd:
                        if medium == 0:
                            os.system('sbatch run_mpipks_sd_sbatch.sh ' + fn_path)
                        elif medium == 1:
                            os.system('sbatch run_mpipks_sd_sbatch_medium.sh ' + fn_path)
                    elif type == FSType.mpipks_mv:
                        if medium == 0:
                            os.system('sbatch run_mpipks_mv_sbatch.sh ' + fn_path)
                        elif medium == 1:
                            os.system('sbatch run_mpipks_mv_sbatch_medium.sh ' + fn_path)
