clear all;

data_path = '../../../../data/cluster/unn/qjx';

sys_id = 0;
task_id = 1;
prop_id = 1;
seed = 0;
mns = 1000000;
num_trajectories = 32;
num_tp_periods = 1000;
num_obs_periods = 1000;
ex_deep = 16;
rk_ns = 10000;

N = 500;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;
dimer_drv_type = 1;
dimer_drv_ampl = 3.4;
dimer_drv_freq = 1;
dimer_drv_phase = 0;
dimer_prm_E = 0;
dimer_prm_U = 0.1125;
dimer_prm_J = 1;
start_type = 0;
start_state = 0;

sys_size = N + 1;

num_runs = 20;


num_bins_x = N+1;
x_min = 0;
x_max = N;
x_shift = (x_max - x_min) / num_bins_x;
x_bins = linspace(x_min + 0.5 * x_shift, x_max - 0.5 * x_shift, num_bins_x);

    
curr_rho_avg = zeros(num_bins_x, 1);

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
        dimer_prm_U, ...
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
        dimer_prm_U, ...
        dimer_prm_J, ...
        start_type, ...
        start_state);
    
    fn = sprintf('%s/mean_evo_%s.txt', path_to_folder, suffix);
    mean_evo = importdata(fn);
    
    x_pdf = zeros(num_bins_x, 1);
    
    for tr_id = 1:num_trajectories
        xs = mean_evo(:, tr_id);
        for p_id = 1 : size(xs, 1)
            x = xs(p_id);
            x_id = floor((x - x_min) / (x_max - x_min + eps) * num_bins_x) + 1;
            x_pdf(x_id) = x_pdf(x_id) + 1;
        end
        
    end
    
    x_pdf = x_pdf / (size(mean_evo, 1) * size(mean_evo, 2) * x_shift);
    
    curr_rho_avg = curr_rho_avg + x_pdf;
end

curr_rho_avg = curr_rho_avg / num_runs;

norm = sum(curr_rho_avg);
norm_diff = 1.0 - norm

states = linspace(1, sys_size, sys_size)';

fig = figure;
hLine = plot(states, curr_rho_avg);
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$|\rho_{n,n}|^2$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);

