clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';

data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 2;
task_id = 7;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 20;
num_tp_periods = 100;
num_obs_periods = 100;
ex_deep = 16;
rk_ns = 10000;

d = 0;
g = 0;

diss_type = 1;
ps_diss_w = 0.05;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
drv_part_1 = 1.0;
drv_part_2 = 1.0;
ps_prm_alpha = 5;
start_type = 0;
start_state = 0;

num_runs = 1;

ampl_begin = 0.05;
ampl_step = 0.05;
ampl_num = 100;
ampls = linspace(ampl_begin, ampl_begin + (ampl_num - 1) * ampl_step, ampl_num);

T_begin = 0.05;
T_step = 0.05;
T_num = 100;
Ts = linspace(T_begin * (drv_part_1 + drv_part_2) * ps_prm_alpha, (T_begin + (T_num - 1) * T_step) * (drv_part_1 + drv_part_2) * ps_prm_alpha, T_num);

lambdas = zeros(ampl_num, T_num);

for T_id = 1:T_num
    
    T = T_begin + (T_id - 1) * T_step
    
    ps_drv_part_1 = drv_part_1 * T;
    ps_drv_part_2 = drv_part_2 * T;
    
    for ampl_id = 1:ampl_num
        
        ampl = ampl_begin + (ampl_id - 1) * ampl_step;
        
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
                ampl, ...
                ps_prm_alpha, ...
                d, ...
                g, ...
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
                ampl, ...
                ps_prm_alpha, ...
                d, ...
                g, ...
                start_type, ...
                start_state);
            
            path = sprintf('%s/lambda_%s.txt', path_to_folder, suffix)
            data = importdata(path);
            
            lambdas(ampl_id, T_id) = lambdas(ampl_id, T_id) + mean(data(num_trajectories / 2 + 1:end));
        end
    end  
end

lambdas = lambdas / num_runs;

fig = figure;
hLine = imagesc(ampls, Ts, lambdas');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$T$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\lambda');
set(gca,'YDir','normal');

suffix = sprintf("s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_var)_prm(%0.4f_%0.4f_%0.4f)", ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    diss_type, ...
    ps_diss_w, ...
    ps_drv_part_1, ...
    ps_drv_part_2, ...
    ps_prm_alpha, ...
    d, ...
    g);

savefig(sprintf('%s/lambda_from_ampl_and_T_%s.fig', home_figures_path, suffix));
