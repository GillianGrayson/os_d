clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';

data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

num_runs = 1;
num_trajectories = 200;

jcs_drv_part_1 = 1.0;
jcs_drv_part_2 = 1.0;
jcs_prm_alpha = 5;

ampl_begin = 0.05;
ampl_step = 0.05;
ampl_num = 100;
ampls = linspace(ampl_begin, ampl_begin + (ampl_num - 1) * ampl_step, ampl_num);

T_begin = 0.05;
T_step = 0.05;
T_num = 100;
Ts = linspace(T_begin * (jcs_drv_part_1 + jcs_drv_part_2) * jcs_prm_alpha, (T_begin + (T_num - 1) * T_step) * (jcs_drv_part_1 + jcs_drv_part_2) * jcs_prm_alpha, T_num);

sys_id = 1;
task_id = 7;
prop_id = 0;

lpn_type = -1;
lpn_delta_s = log10(1.0e-3);
lpn_delta_f_high = log10(1.0e-1);
lpn_delta_f_low = log10(1.0e-5);


ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 100;
num_obs_periods = 100;

num_random_obs = 1;
random_obs_seed = 100;
random_obs_type = 2;

N = 300;

diss_type = 0;

start_type = 0;
start_state = 0;

seed = 0;
mns = 1000000;

lambdas = zeros(ampl_num, T_num);

for T_id = 1:T_num
    
    T = T_begin + T_step * (T_id - 1)
    
    for ampl_id = 1:ampl_num
        
        ampl = ampl_begin + ampl_step * (ampl_id - 1);
        
        for run_id = 1:num_runs
            
            ss = (run_id - 1) * num_trajectories;
            
            path_to_folder = sprintf('%s/main_%d_%d_%d/lpn_%d_%0.4f_%0.4f_%0.4f/run_%d_%d_%d_%d/obs_%d_%d_%d/N_%d/diss_%d/drv_%0.4f_%0.4f_%0.4f/prm_%0.4f/start_%d_%d/ss_%d', ...
                data_path, ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                lpn_type, ...
                lpn_delta_s, ...
                lpn_delta_f_high, ...
                lpn_delta_f_low, ...
                ex_deep, ...
                rk_ns, ...
                num_tp_periods, ...
                num_obs_periods, ...
                num_random_obs, ...
                random_obs_seed, ...
                random_obs_type, ...
                N, ...
                diss_type, ...
                jcs_drv_part_1 * T, ...
                jcs_drv_part_2 * T, ...
                ampl, ...
                jcs_prm_alpha, ...
                start_type, ...
                start_state, ...
                ss);
            
            suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                ss, ...
                mns, ...
                N, ...
                diss_type, ...
                0.1, ...
                0.0, ...
                jcs_drv_part_1 * T, ...
                jcs_drv_part_2 * T, ...
                ampl, ...
                jcs_prm_alpha, ...
                start_type, ...
                start_state, ...
                lpn_type, ...
                lpn_delta_s, ...
                lpn_delta_f_high, ...
                lpn_delta_f_low);
            
            
            path = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
            data = importdata(path);
            
            lambdas(ampl_id, T_id) = lambdas(ampl_id, T_id) + mean(data(num_trajectories / 2 + 1:end));
        end
    end
end

lambdas = lambdas / num_runs;

suffix = sprintf('lpn(%d_%0.4f_%0.4f_%0.4f)_N(%d)_drv(var_var_var)_prm(%0.4f)_start(%d_%d)_time(%d_%d)_numTraj(%d)', ...
            lpn_type, ...
            lpn_delta_s, ...
            lpn_delta_f_high, ...
            lpn_delta_f_low,...
			N, ...
            jcs_prm_alpha, ...
            start_type, ...
            start_state, ...
			num_tp_periods, ...
            num_obs_periods, ...
			num_trajectories/2);

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

savefig(sprintf('%s/lambda_from_ampl_and_T_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_from_ampl_and_T_%s.pdf', home_figures_path, suffix));
