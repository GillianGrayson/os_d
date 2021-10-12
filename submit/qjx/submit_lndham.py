import pathlib
from Infrastructure.file_system import *
import os.path
import numpy as np

type = FSType.mpipks_sd

num_runs = 1

medium = 1

seed_start = 1
seed_shift = 1
seed_num = 100

for seed_id in range(0, seed_num):
    lndham_seed = seed_start + seed_id * seed_shift

    print('mbl_seed: ' + str(lndham_seed))

    sys_id = 4
    task_id = 7  # 7
    prop_id = 0
    is_debug = 0
    is_pp = 0
    init_fn = ''
    path = ''
    seed = 0
    mns = 1000000
    num_threads = 1
    num_trajectories = 200  # 200
    num_tp_periods = 500
    num_obs_periods = 500  # 500
    ex_deep = 16
    rk_ns = 10000

    num_random_obs = 1
    random_obs_seed = lndham_seed
    random_obs_mns = 1000000
    random_obs_type = 2
    lpn_type = -1
    lpn_eps_deep = 100
    lpn_eps_error = 1.0e-10
    lpn_eps_high = 10
    lpn_eps_low = -10
    lpn_delta_s = 1.0e-6
    lpn_delta_f_high = 1.0e-6
    lpn_delta_f_low = 1.0e-6
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
    lndham_N = 100
    lndham_seed = lndham_seed
    lndham_num_seeds = 1000000
    lndham_alpha = 0.5
    lndham_T = 1.0
    start_type = 0
    start_state = 0
    cd_dim = 1
    cd_eps = 1.0e-8
    deep_num_steps = 1000
    jump = 0  # 0
    jumps_counts = 0  # 0

    diss_phase_str = str(format(diss_phase, '0.4f'))
    diss_gamma_str = str(format(diss_gamma, '0.4f'))

    lndham_alpha_str = str(format(lndham_alpha, '0.4f'))
    lndham_T_str = str(format(lndham_T, '0.4f'))

    lpn_delta_s_str = str(format(np.log10(lpn_delta_s), '0.4f'))
    lpn_delta_f_high_str = str(format(np.log10(lpn_delta_f_high), '0.4f'))
    lpn_delta_f_low_str = str(format(np.log10(lpn_delta_f_low), '0.4f'))

    start_seed = 0
    finish_seed = num_runs * num_trajectories
    step_seed = num_trajectories

    print('start_seed: ' + str(start_seed))
    print('finish_seed: ' + str(finish_seed))
    print('step_seed: ' + str(step_seed))

    local_path = f"/main_{sys_id}_{task_id}_{prop_id}"

    if task_id == 7:
        local_path += f"/lpn_{lpn_type}_{lpn_delta_s_str}"

    local_path += f"/run_{ex_deep}_{rk_ns}_{num_tp_periods}_{num_obs_periods}" + \
                  f"/N_{lndham_N}_alpha_{lndham_alpha_str}" + \
                  f'/seed_{lndham_seed}'

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
        file_params.write('lndham_N ' + str(lndham_N) + '\n')
        file_params.write('lndham_seed ' + str(lndham_seed) + '\n')
        file_params.write('lndham_num_seeds ' + str(lndham_num_seeds) + '\n')
        file_params.write('lndham_alpha ' + str(lndham_alpha) + '\n')
        file_params.write('lndham_T ' + str(lndham_T) + '\n')
        file_params.write('start_type ' + str(start_type) + '\n')
        file_params.write('start_state ' + str(start_state) + '\n')
        file_params.write('cd_dim ' + str(cd_dim) + '\n')
        file_params.write('cd_eps ' + str(cd_eps) + '\n')
        file_params.write('deep_num_steps ' + str(deep_num_steps) + '\n')
        file_params.write('jump ' + str(jump) + '\n')
        file_params.write('jumps_counts ' + str(jumps_counts))

        file_params.close()

        fn_suffix = \
            'setup(' + str(sys_id) + '_' + str(task_id) + '_' + str(prop_id) + ')_' + \
            'rnd(' + str(ss) + '_' + str(mns) + ')_' + \
            'N(' + str(lndham_N) + ')_' + \
            'rnd(' + str(lndham_seed) + ')_' + \
            'alpha(' + lndham_alpha_str + ')_' + \
            'T(' + lndham_T_str  + ')'

        if task_id == 7:
            fn_suffix += '_lpn(' + str(
                lpn_type) + '_' + lpn_delta_s_str + '_' + lpn_delta_f_high_str + '_' + lpn_delta_f_low_str + ')'

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
