clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'qjx_results';

data_path = sprintf('%s/%s', data_path, prefix);

sys_id = 0;
task_id = 0;
prop_id = 1;
is_debug = 0;
is_pp = 1;
init_fn = '';
path = '';
num_threads = 1;
qj_deep = 16;
qj_num_tp_periods = 1000;
qj_num_obs_periods = 1000;
qj_num_trajectories = 2;
qj_seed = 0;
qj_mns = 1000000;

type_lpn = 0;
eps_lpn = 1.0e-3;
delta_up_lpn = 1.0e-2;
delta_down_lpn = 1.0e-13;
is_obs_dump = 1;
is_adr_dump_sep = 0;
is_adr_dump_avg = 0;
is_evo_dump_sep = 1;
is_evo_dump_avg = 0;
dump_type = 0;
num_dumps = 1;
N = 100;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;
start_type = 0;
start_state = 0;
prm_E = 0.0;
prm_U = 0;
prm_J = 1.0;

sys_size = N + 1;

num_runs = 100;

U_begin = 0.01;
U_step = 0.01;
U_num = 75;

for dump_id = 101 : 100 : 1001

	Us = zeros(U_num, 1);
	lambdas = zeros(U_num, 1);

	for U_id = 1:U_num
		
		U_id = U_id
		U = U_begin + U_step * (U_id - 1);
		Us(U_id) = U;
	   
		curr_lambda_avg = 0;
	   
		for run_id = 1:num_runs
			
			ss = (run_id - 1) * qj_num_trajectories;
			
			path_to_folder = sprintf('%s/main_%d_%d_%d/qj_%d_%d_%d/N_%d/diss_%d_%0.4f_%0.4f/drv_%d_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
				data_path, ...
				sys_id, ...
				prop_id, ...
				task_id, ...
				qj_deep, ...
				qj_num_tp_periods, ...
				qj_num_obs_periods, ...
				N, ...
				diss_type, ...
				diss_gamma, ...
				diss_phase, ...
				drv_type, ...
				drv_ampl, ...
				drv_freq, ...
				drv_phase, ...
				prm_E, ...
				U, ...
				prm_J, ...
				start_type, ...
				start_state, ...
				ss);
			
			suffix = sprintf('rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
				ss, ...
				qj_mns, ...
				N, ...
				diss_type, ...
				diss_gamma, ...
				diss_phase, ...
				drv_type, ...
				drv_ampl, ...
				drv_freq, ...
				drv_phase, ...
				prm_E, ...
				U, ...
				prm_J, ...
				start_type, ...
				start_state);
		   
		  
			path = sprintf('%s/lambda_evo_%s.txt', path_to_folder, suffix);
			data = importdata(path);
		   
			curr_lambda_avg = curr_lambda_avg + data(dump_id, 2);
		end
	   
		curr_lambda_avg = curr_lambda_avg / num_runs;
	   
		lambdas(U_id) = curr_lambda_avg;
		
	end

	fig = figure;
	hLine = plot(Us, lambdas);
	set(gca, 'FontSize', 30);
	xlabel('$U$', 'Interpreter', 'latex');
	set(gca, 'FontSize', 30);
	ylabel('$\lambda$', 'Interpreter', 'latex');

	suffix = sprintf('qj_dump(%d)_ns(%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
				dump_id, ...
				num_runs * qj_num_trajectories, ...
				N, ...
				diss_type, ...
				diss_gamma, ...
				diss_phase, ...
				drv_type, ...
				drv_ampl, ...
				drv_freq, ...
				drv_phase, ...
				prm_E, ...
				prm_U, ...
				prm_J, ...
				start_type, ...
				start_state);
			
	savefig(sprintf('%s/lambda_%s.fig', home_figures_path, suffix));

end