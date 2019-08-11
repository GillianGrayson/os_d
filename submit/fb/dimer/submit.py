import pathlib
from Infrastructure.file_system import *
import os.path

type = FSType.mpipks_sd

U_start = 0.005
U_shift = 0.005
U_num = 150

ampl_start = 3.4
ampl_shift = 0.05
ampl_num = 1

for U_id in range(0, U_num):

    U = U_start + U_id * U_shift

    for ampl_id in range(0, ampl_num):
        ampl = ampl_start + ampl_id * ampl_shift

        print('U: ' + str(U))
        print('A: ' + str(ampl))

        task = 2
        debug = 0
        issmtx = 1
        ipp = 1
        N = 100
        diss_type = 0
        diss_gamma = 0.1
        diss_phase = 0.0
        drv_type = 1
        drv_ampl = 3.4
        drv_freq = 1.0
        drv_phase = 0.0
        prm_E = 0.0
        prm_U = U
        prm_J = 1.0
        num_steps = 10000
        num_periods_trans = 100
        num_periods_obser = 0
        seed = 1
        max_num_seeds = 1000000
        int_ist = 1
        int_isi = 0
        int_dt = 0
        int_dn = 1

        diss_gamma_str = str(format(diss_gamma, '0.4f'))
        diss_phase_str = str(format(diss_phase, '0.4f'))

        drv_ampl_str = str(format(drv_ampl, '0.4f'))
        drv_freq_str = str(format(drv_freq, '0.4f'))
        drv_phase_str = str(format(drv_phase, '0.4f'))

        prm_E_str = str(format(prm_E, '0.4f'))
        prm_U_str = str(format(prm_U, '0.4f'))
        prm_J_str = str(format(prm_J, '0.4f'))

        local_path = \
            '/task_' + str(task) + \
            '/int_' + str(num_steps) + '_' + str(num_periods_trans) + '_' + str(num_periods_obser) + \
            '/N_' + str(N) + \
            '/diss_' + str(diss_type) + '_' + diss_gamma_str + '_' + diss_phase_str + \
            '/drv_' + str(drv_type) + '_' + drv_ampl_str + '_' + drv_freq_str + '_' + drv_phase_str + \
            '/prm_' + prm_E_str + '_' + prm_U_str + '_' + prm_J_str + \
			'/seed_' + str(seed)

        fn_path = get_dir(local_path, type)
        pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

        file_config = open(fn_path + '/config.txt', 'w')
        file_config.write('task ' + str(task) + '\n')
        file_config.write('debug ' + str(debug) + '\n')
        file_config.write('issmtx ' + str(issmtx) + '\n')
        file_config.write('ipp ' + str(ipp) + '\n')
        file_config.write('N ' + str(N) + '\n')
        file_config.write('diss_type ' + str(diss_type) + '\n')
        file_config.write('diss_gamma ' + str(diss_gamma) + '\n')
        file_config.write('diss_phase ' + str(diss_phase) + '\n')
        file_config.write('drv_type ' + str(drv_type) + '\n')
        file_config.write('drv_ampl ' + str(drv_ampl) + '\n')
        file_config.write('drv_freq ' + str(drv_freq) + '\n')
        file_config.write('drv_phase ' + str(drv_phase) + '\n')
        file_config.write('prm_E ' + str(prm_E) + '\n')
        file_config.write('prm_U ' + str(prm_U) + '\n')
        file_config.write('prm_J ' + str(prm_J) + '\n')
        file_config.write('num_steps ' + str(num_steps) + '\n')
        file_config.write('num_periods_trans ' + str(num_periods_trans) + '\n')
        file_config.write('num_periods_obser ' + str(num_periods_obser) + '\n')
        file_config.write('seed ' + str(seed) + '\n')
        file_config.write('max_num_seeds ' + str(max_num_seeds) + '\n')
        file_config.write('int_ist ' + str(int_ist) + '\n')
        file_config.write('int_isi ' + str(int_isi) + '\n')
        file_config.write('int_dt ' + str(int_dt) + '\n')
        file_config.write('int_dn ' + str(int_dn) + '\n')
        file_config.close()

        if task == 0:
            fn_test = fn_path + '/rho.txt'
        elif task == 1 or task == 2:
            fn_test = fn_path + '/floquet_evals.txt'

        if not os.path.isfile(fn_test):
            if type == FSType.unn:
                os.system('sbatch run_unn.sh ' + fn_path)
            elif type == FSType.mpipks_sd:
                os.system('sbatch run_mpipks.sh ' + fn_path)
