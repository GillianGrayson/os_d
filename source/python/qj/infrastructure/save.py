from infrastructure.path import *


def create_file(config):
    path = get_data_path(config)

    f = open(path + '/config.txt', 'w')

    f.write('sys_id ' + str(config.run.quantum_system.value) + '\n')
    f.write('task_id ' + str(config.run.process.value) + '\n')
    f.write('prop_id ' + str(config.run.propagation.value) + '\n')
    f.write('ex_deep ' + str(config.details.step_metrics) + '\n')
    f.write('rk_ns ' + str(config.details.step_metrics) + '\n')
    f.write('num_tp_periods ' + str(config.details.num_periods_trans) + '\n')
    f.write('num_obs_periods ' + str(config.details.num_periods_obser) + '\n')
    f.write('seed ' + str(config.random.seed) + '\n')
    f.write('mns ' + str(config.random.num_seeds) + '\n')
    f.write('is_debug ' + str(config.auxiliary['is_debug']) + '\n')
    f.write('is_pp ' + str(config.auxiliary['is_pp']) + '\n')
    f.write('init_fn ' + str(config.auxiliary['init_fn']) + '\n')
    f.write('path ' + str(config.auxiliary['path']) + '\n')
    f.write('num_threads ' + str(config.auxiliary['num_threads']) + '\n')
    f.write('num_trajectories ' + str(config.auxiliary['num_trajectories']) + '\n')

    f.close()

    f = open(path + '/params.txt', 'w')

    f.write('lpn_type ' + str(config.auxiliary['lpn_type']) + '\n')
    f.write('lpn_eps_deep ' + str(config.auxiliary['lpn_eps_deep']) + '\n')
    f.write('lpn_eps_high ' + str(config.auxiliary['lpn_eps_high']) + '\n')
    f.write('lpn_eps_low ' + str(config.auxiliary['lpn_eps_low']) + '\n')
    f.write('lpn_delta_s_high ' + str(config.auxiliary['lpn_delta_s_high']) + '\n')
    f.write('lpn_delta_s_low ' + str(config.auxiliary['lpn_delta_s_low']) + '\n')
    f.write('lpn_delta_f_high ' + str(config.auxiliary['lpn_delta_f_high']) + '\n')
    f.write('lpn_delta_f_low ' + str(config.auxiliary['lpn_delta_f_low']) + '\n')
    f.write('dump_obs ' + str(config.auxiliary['dump_obs']) + '\n')
    f.write('dump_adr_sep ' + str(config.auxiliary['dump_adr_sep']) + '\n')
    f.write('dump_adr_avg ' + str(config.auxiliary['dump_adr_avg']) + '\n')
    f.write('dump_evo_sep ' + str(config.auxiliary['dump_evo_sep']) + '\n')
    f.write('dump_evo_avg ' + str(config.auxiliary['dump_evo_avg']) + '\n')
    f.write('dump_type ' + str(config.auxiliary['dump_type']) + '\n')
    f.write('dump_num ' + str(config.auxiliary['dump_num']) + '\n')
    f.write('cd_dim ' + str(config.auxiliary['cd_dim']) + '\n')
    f.write('cd_eps ' + str(config.auxiliary['cd_eps']) + '\n')
    f.write('deep_num_steps ' + str(config.auxiliary['deep_num_steps']) + '\n')
    f.write('jump ' + str(config.auxiliary['jump']) + '\n')
    f.write('jumps_counts ' + str(config.auxiliary['jumps_counts']) + '\n')

    params_size = config.params.size
    for key, value in params_size.items():
        f.write(key + ' ' + str(value) + '\n')
    params_dissipation = config.params.dissipation
    for key, value in params_dissipation.items():
        f.write(key + ' ' + str(value) + '\n')
    params_driving = config.params.driving
    for key, value in params_driving.items():
        f.write(key + ' ' + str(value) + '\n')
    params_model = config.params.model
    for key, value in params_model.items():
        f.write(key + ' ' + str(value) + '\n')
    params_init_cond = config.params.init_cond
    for key, value in params_init_cond.items():
        f.write(key + ' ' + str(value) + '\n')

    f.close()