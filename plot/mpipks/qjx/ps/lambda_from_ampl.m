clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';

data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';


sys_id = 2; 
task_id = 0; 
prop_id = 0; 
seed = 0; 
mns = 1000000;
num_trajectories = 100;
num_tp_periods = 100;
num_obs_periods = 100; 
ex_deep = 16;
rk_ns = 10000;

diss_type = 1;
ps_diss_w = 0.05;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
ps_drv_part_1 = 0.98;
ps_drv_part_2 = 1.00;
ps_drv_ampl = 3.2;
ps_prm_alpha = 5;
ps_prm_d = 0.;
ps_prm_g = 0.;
start_type = 0;
start_state = 0;

num_runs = 10;

num_points = 100;
lambdas = zeros(num_points, 1);
params = zeros(num_points, 1);
for param_id = 1:num_points
    
    ps_drv_ampl = 0.05 + (param_id - 1) * 0.05
    params(param_id) = ps_drv_ampl;
    
    for run_id = 1:num_runs
        
        ss = (run_id - 1) * num_trajectories;
        
        path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/N_%d_%d/diss_%d_%0.4f/drv_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
            data_path, ...
            sys_id, ...
            task_id, ...
            prop_id, ...
            ex_deep, ...
            rk_ns, ...
            num_tp_periods, ...
            num_obs_periods, ...
            ps_num_spins, ...
            ps_num_photons_states, ...
            diss_type, ...
            ps_diss_w, ...
            ps_drv_part_1, ...
            ps_drv_part_2, ...
            ps_drv_ampl, ...
            ps_prm_alpha, ...
            ps_prm_d, ...
            ps_prm_g, ...
            start_type, ...
            start_state, ...
            ss);
        
        suffix = sprintf("rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
            ss, ...
            mns, ...
            ps_num_spins, ...
            ps_num_photons_states, ...
            diss_type, ...
            ps_diss_w, ...
            ps_drv_part_1, ...
            ps_drv_part_2, ...
            ps_drv_ampl, ...
            ps_prm_alpha, ...
            ps_prm_d, ...
            ps_prm_g, ...
            start_type, ...
            start_state);
        
        path = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
        data = importdata(path);
        
        lambdas(param_id) = lambdas(param_id) + mean(data(2:end));
    end
    
    lambdas(param_id) = lambdas(param_id) / num_runs;
    
    
end

fig = figure;
hLine = plot(params, lambdas);
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');
hold all;

suffix = sprintf("s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_var)_prm(%0.4f_%0.4f_%0.4f)", ...
	ps_num_spins, ...
	ps_num_photons_states, ...
	diss_type, ...
	ps_diss_w, ...
	ps_drv_part_1, ...
	ps_drv_part_2, ...
	ps_prm_alpha, ...
	ps_prm_d, ...
	ps_prm_g);

savefig(sprintf('%s/lambda_from_ampl_%s.fig', home_figures_path, suffix));
