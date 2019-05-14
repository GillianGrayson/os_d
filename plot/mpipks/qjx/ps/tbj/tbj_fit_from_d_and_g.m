clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

diss = 0;

sys_id = 2;
task_id = 1;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 1;
num_tp_periods = 100;
num_obs_periods = 20000;
ex_deep = 16;
rk_ns = 10000;

diss_type = 1;
ps_diss_w = 0.05;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
ps_drv_part_1 = 0.98;
ps_drv_part_2 = 1.00;
ps_drv_ampl = 2.5;
ps_prm_alpha = 5;
start_type = 0;
start_state = 0;

num_runs = 100;

param_start = 0.0;
param_shift = 0.1;
param_num = 100;

params = zeros(param_num, 1);
alphas = zeros(param_num, 1);
decades = zeros(param_num, 1);

for param_id = 1:param_num
    
    param = param_start + param_shift * param_id
    params(param_id) = param;
    
    ps_prm_d = param;
    ps_prm_g = param;
    
    bin_begin = 1e-10;
    num_decades = 15;
    bin_end = bin_begin * 10.^num_decades;
    num_bin_per_decade = 10;
    num_bins = num_bin_per_decade * num_decades;
    bin_borders = zeros(num_bins + 1, 1);
    bin_centers = zeros(num_bins, 1);
    for bin_id = 1 : num_bins + 1
        bin_borders(bin_id) = bin_begin * 10.^((bin_id - 1) / num_bin_per_decade);
        if (bin_id <= num_bins)
            bin_centers(bin_id) = bin_begin * 10.^((bin_id - 1 + 0.5) / num_bin_per_decade);
        end
    end
    bin_diff = diff(bin_borders);
    
    non_inc_count = 0;
    curr_pdf = zeros(num_bins, 1);
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
        
        for trajectory_id = 1:num_trajectories
            
            path = sprintf('%s/jump_times_%d_%s.txt', path_to_folder, trajectory_id - 1, suffix);
            data = importdata(path);
            path = sprintf('%s/diss_types_%d_%s.txt', path_to_folder, trajectory_id - 1, suffix);
            diss_types = importdata(path);
            if diss ~= -1
                indexes = find(diss_types ~= diss);
                num_diss_before = size(data, 1);
                data(indexes) = [];
                num_diss_after = size(data, 1);
            end
            curr_diff = diff(data);
            curr_size = size(curr_diff, 1);
            
            for jump_diff_id = 1:curr_size
                
                curr_tmp = curr_diff(jump_diff_id);
                
                if ((curr_tmp >= bin_begin) && (curr_tmp <= bin_end))
                    bin_id = floor((log10(curr_tmp) - log10(bin_begin)) * num_bins / (log10(bin_end) - log10(bin_begin) + eps)) + 1;
                    curr_pdf(bin_id) = curr_pdf(bin_id) + 1;
                else
                    non_inc_count = non_inc_count + 1;
                end
                
            end
        end
        
    end
    
    non_inc_count = non_inc_count
    norm = sum(curr_pdf);
    
    for bin_id = 1 : num_bins
        curr_pdf(bin_id) = curr_pdf(bin_id) / (norm * bin_diff(bin_id));
    end
    
    norm_check = 0;
    for bin_id = 1 : num_bins
        norm_check = norm_check + curr_pdf(bin_id) * bin_diff(bin_id);
    end
    norm_check = norm_check
    
    [x_min, x_max] = find_range(bin_centers, curr_pdf);
    [alpha, coef, R2, yy, R2_, cnt] = powerlaw_regression(bin_centers, curr_pdf, x_min, x_max);
    left = find(bin_centers==x_min);
    right = find(bin_centers==x_max);
    
    if ~isnan(R2)
        alphas(param_id) = abs(alpha);
        decades(param_id) = log10(x_max) - log10(x_min);
    else
        alphas(param_id) = 0;
        decades(param_id) = 0;
    end
   
end

suffix_save = sprintf("diss(%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_var_var)", ...
    diss, ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    diss_type, ...
    ps_diss_w, ...
    ps_drv_part_1, ...
    ps_drv_part_2, ...
    ps_drv_ampl, ...
    ps_prm_alpha);

fig = figure;
hLine = plot(params, alphas, 'LineWidth', 3);
set(gca, 'FontSize', 30);
xlabel('$d=g$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\alpha$', 'Interpreter', 'latex');

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/tbj_alphas_from_d_and_g_%s.pdf', home_figures_path, suffix_save));
savefig(gcf, sprintf('%s/tbj_alphas_from_d_and_g_%s.fig', home_figures_path, suffix_save));
close(gcf)

fig = figure;
hLine = plot(params, decades, 'LineWidth', 3);
set(gca, 'FontSize', 30);
xlabel('$d=g$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$decades$', 'Interpreter', 'latex');

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/tbj_decades_from_d_and_g_%s.pdf', home_figures_path, suffix_save));
savefig(gcf, sprintf('%s/tbj_decades_from_d_and_g_%s.fig', home_figures_path, suffix_save));
close(gcf)