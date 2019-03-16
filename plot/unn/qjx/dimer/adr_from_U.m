clear all;

data_path = '../../../../data/cluster/unn/qjx';

sys_id = 0;
task_id = 1;
prop_id = 1;
seed = 0;
mns = 1000000;
num_trajectories = 32;
num_tp_periods = 1000;
num_obs_periods = 10000;
ex_deep = 16;
rk_ns = 10000;

N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;
dimer_drv_type = 1; 
dimer_drv_ampl = 4.2; 
dimer_drv_freq = 1; 
dimer_drv_phase = 0; 
dimer_prm_E = 0;
dimer_prm_U = 0.02;
dimer_prm_J = 1; 
start_type = 0;
start_state = 0;

sys_size = N + 1;

num_runs = 1;

U_begin = 0.01;
U_step = 0.01;
U_num = 100;

Us = zeros(U_num, 1);
states = linspace(1, sys_size, sys_size);

rho = zeros(U_num, sys_size);

for U_id = 1:U_num
	
	U_id = U_id
	U = U_begin + U_step * (U_id - 1);
	Us(U_id) = U;
   
	curr_rho_avg = zeros(sys_size,1);
   
	for run_id = 1:num_runs
        
        ss = (run_id - 1) * num_trajectories;
        
		path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/N_%d/diss_%d_%0.4f_%0.4f/drv_%d_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
			data_path, ...
            sys_id, ...
			task_id, ...
			prop_id, ...
            ex_deep, ...
			rk_ns, ...
            num_tp_periods, ...
            num_obs_periods, ...
            N, ...
            diss_type, ...
            diss_gamma, ...
            diss_phase, ...
            dimer_drv_type, ...
            dimer_drv_ampl, ...
            dimer_drv_freq, ...
            dimer_drv_phase, ...
            dimer_prm_E, ...
            U, ...
            dimer_prm_J, ...
            start_type, ...
            start_state, ...
            ss);
        
    suffix = sprintf('rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
            ss, ...
            mns, ...
            N, ...
            diss_type, ...
            diss_gamma, ...
            diss_phase, ...
            dimer_drv_type, ...
            dimer_drv_ampl, ...
            dimer_drv_freq, ...
            dimer_drv_phase, ...
            dimer_prm_E, ...
            U, ...
            dimer_prm_J, ...
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

target_U_id = 2;
target_U = Us(target_U_id)
target_rho = rho(:, target_U_id);

fig = figure;
hLine = plot(states, target_rho);
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$PDF$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);

