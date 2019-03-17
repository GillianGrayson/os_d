clear all;

data_path = '../../../../data/cluster/unn/qjx';

eps = 1.0e-8;

adaptive_axes = 0;

num_runs = 20;

sys_id = 0;
task_id = 1;
prop_id = 1;
seed = 0;
mns = 1000000;
num_trajectories = 32;
num_tp_periods = 100;
num_obs_periods = 100;
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

rho_total = zeros(sys_size);
count_rho = 0;

for run_id = 1:num_runs
    
    ss = 0 + (run_id - 1) * num_trajectories;
    
    for tr_id = 1:num_trajectories
        
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
        
        suffix = sprintf('%d_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
            tr_id - 1, ...
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
        
        fn = sprintf('%s/phi_evo_%s.txt', path_to_folder, suffix);
        phi_evo_raw = importdata(fn);
        phi_evo = complex(phi_evo_raw(:, 1), phi_evo_raw(:, 2));
        
        curr_phi = phi_evo(1:sys_size);
        curr_norm = sum(abs(curr_phi).^2);
        curr_rho = zeros(sys_size);
        for i = 1:sys_size
            for j = 1:sys_size
                curr_rho(i,j) = dot(curr_phi(i), curr_phi(j)) / curr_norm;
            end
        end
        norm_check = 1 - trace(curr_rho)
        count_rho = count_rho + 1
        
        rho_total = rho_total + curr_rho;
        
    end
end

rho_total = rho_total / count_rho;
norm_check_result = 1 - trace(curr_rho)

fig = figure;
propertyeditor(fig);

states = linspace(1, sys_size, sys_size)';
hLine = imagesc(states, states, abs(rho_total));
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$m$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
h.Label.Interpreter = 'latex';
title(h, '$|\rho_{n,m}|$', 'FontSize', 33, 'Interpreter', 'latex');
set(gca,'YDir','normal');
hold all;