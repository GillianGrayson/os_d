import pathlib
from Infrastructure.file_system import *
import os.path
import math
import numpy as np

type = FSType.unn

drv_ampl_start = 4.0
drv_ampl_shift = 0.1
drv_ampl_num = 1

T_start = 0.1
T_shift = 0.1
T_num = 1

for T_id in range(0, T_num):

	T = T_start + T_id * T_shift

	for drv_ampl_id in range(0, drv_ampl_num):

		curr_drv_ampl = drv_ampl_start + drv_ampl_id * drv_ampl_shift

		print('curr_drv_ampl: ' + str(curr_drv_ampl))

		task = 1
		debug = 0
		issmtx = 1
		ipp = 1
		init_file = ''
		path = ''
		N = 200
		dt = 1
		dp = 0.0
		g = 0.1
		g_add = 0.1
		drv_ampl = curr_drv_ampl
		prm_alpha = 5.0
		num_steps_t_0 = 10000
		num_steps_t_1 = 10000
		num_periods_trans = 100
		num_periods_obser = 1
		t_0 = 0.98 * T
		t_1 = 1.00 * T
		seed = 1
		max_num_seeds = 1000000
		int_ist = 1
		int_isi = 0
		int_dt = 0
		int_dn = 100

		drv_ampl_str = str(format(drv_ampl, '0.4f'))
		prm_alpha_str = str(format(prm_alpha, '0.4f'))

		dp_str = str(format(dp, '0.4f'))
		g_str = str(format(g, '0.4f'))
		g_add_str = str(format(g_add, '0.4f'))

		t_0_str = str(format(t_0, '0.4f'))
		t_1_str = str(format(t_1, '0.4f'))


		local_path = \
			'/main_' + str(task) + \
			'/N_' + str(N) + \
			'/prm_' + drv_ampl_str + '_' + prm_alpha_str + '_' + t_0_str + '_' + t_1_str + \
			'/diss_' + str(dt) + '_' + dp_str + '_' + g_str + '_' + g_add_str

		fn_path = get_dir(local_path, type)
		pathlib.Path(fn_path).mkdir(parents=True, exist_ok=True)

		file_config = open(fn_path + '/config.txt', 'w')
		file_config.write('task ' + str(task) + '\n')
		file_config.write('debug ' + str(debug) + '\n')
		file_config.write('issmtx ' + str(issmtx) + '\n')
		file_config.write('ipp ' + str(ipp) + '\n')
		file_config.write('init_file ' + str(init_file) + '\n')
		file_config.write('path ' + str(path) + '\n')
		file_config.write('N ' + str(N) + '\n')
		file_config.write('dt ' + str(dt) + '\n')
		file_config.write('dp ' + str(dp) + '\n')
		file_config.write('g ' + str(g) + '\n')
		file_config.write('g_add ' + str(g_add) + '\n')
		file_config.write('drv_ampl ' + str(drv_ampl) + '\n')
		file_config.write('prm_alpha ' + str(prm_alpha) + '\n')
		file_config.write('num_steps_t_0 ' + str(num_steps_t_0) + '\n')
		file_config.write('num_steps_t_1 ' + str(num_steps_t_1) + '\n')
		file_config.write('num_periods_trans ' + str(num_periods_trans) + '\n')
		file_config.write('num_periods_obser ' + str(num_periods_obser) + '\n')
		file_config.write('t_0 ' + str(t_0) + '\n')
		file_config.write('t_1 ' + str(t_1) + '\n')
		file_config.write('seed ' + str(seed) + '\n')
		file_config.write('max_num_seeds ' + str(max_num_seeds) + '\n')
		file_config.write('int_ist ' + str(int_ist) + '\n')
		file_config.write('int_isi ' + str(int_isi) + '\n')
		file_config.write('int_dt ' + str(int_dt) + '\n')
		file_config.write('int_dn ' + str(int_dn) + '\n')
		file_config.close()

		fn_test = fn_path + '/rho.txt'

		if not os.path.isfile(fn_test):
			if type == FSType.unn:
				os.system('sbatch run_unn.sh ' + fn_path)
			elif type == FSType.mpipks:
				os.system('qsub run_mpipks.sh ' + fn_path)
