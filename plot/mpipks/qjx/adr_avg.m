clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'qjx_results';

data_path = sprintf('%s/%s', data_path, prefix);

sys_id = 0;
task_id = 1;
is_debug = 0;
is_pp = 1;
init_fn = '';
path = '';
num_threads = 1;
qj_deep = 16;
qj_num_tp_periods = 100;
qj_num_obs_periods = 1;
qj_num_trajectories = 1000;
qj_seed = 0;
qj_mns = 1000000;

type_lpn = 0;
eps_lpn = 1.0e-3;
delta_up_lpn = 1.0e-2;
delta_down_lpn = 1.0e-13;
is_obs_dump = 1;
is_adr_dump_sep = 0;
is_adr_dump_avg = 1;
is_evo_dump_sep = 0;
is_evo_dump_avg = 0;
dump_type = 0;
num_dumps = 1;
N = 100;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;
start_type = 0;
start_state = 0;
prm_E = 1.0;
prm_U = 0;
prm_J = 1.0;

sys_size = N + 1;

num_runs = 100;

U_begin = 0.01;
U_step = 0.01;
U_num = 75;

Us = zeros(U_num, 1);
states = linspace(1, sys_size, sys_size);

rho = zeros(U_num, sys_size);

for U_id = 1:U_num
	
	U_id = U_id
	U = U_begin + U_step * (U_id - 1);
	Us(U_id) = U;
   
	curr_rho_avg = zeros(sys_size,1);
   
	for run_id = 1:num_runs
        
        ss = (run_id - 1) * qj_num_trajectories;
        
		path_to_folder = sprintf('%s/main_%d_%d/qj_%d_%d_%d/N_%d/diss_%d_%0.4f_%0.4f/drv_%d_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
			data_path, ...
            sys_id, ...
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
        
        suffix = sprintf('qjrnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
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
       
      
		path = sprintf('%s/adr_avg_%s.txt', path_to_folder, suffix);
		data = importdata(path);
       
		curr_rho_avg = curr_rho_avg + data;
	end
   
	curr_rho_avg = curr_rho_avg / num_runs;
   
	norm = sum(curr_rho_avg);
    norm_diff = 1.0 - norm
   
	rho(U_id, :) = curr_rho_avg;
    
end

for U_id = 1:U_num
	curr_max = max(rho(U_id, :));
	rho(U_id, :) = rho(U_id, :) / curr_max; 
end

rho = rho';

fig = figure;
hLine = imagesc(Us, states, rho + eps);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\rho_{n,n}');
set(gca,'YDir','normal');

suffix = sprintf('qj(%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', ...
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
		
savefig(sprintf('%s/adr_avg_%s.fig', home_figures_path, suffix));