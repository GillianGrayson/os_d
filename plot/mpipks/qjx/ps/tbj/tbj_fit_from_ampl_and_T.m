clear all;

addpath('/home/ivanchen/yusipov/os_lnd/source/matlab/lib')

home_figures_path = '/home/ivanchen/yusipov/os_d/figures';

data_path = '/data/condmat/ivanchen/yusipov/os_d/qjx';

diss = 0;

sys_id = 2;
task_id = 1;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 20;
num_tp_periods = 10;
num_obs_periods = 20000;
ex_deep = 16;
rk_ns = 10000;

num_random_obs = 1;
random_obs_seed = 100;
random_obs_type = 2;


diss_type = 1;
ps_diss_w = 0.05;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 300;
ps_drv_part_1 = 1.00;
ps_drv_part_2 = 1.00;
ps_prm_alpha = 5;
start_type = 0;
start_state = 0;

ampl_start = 0.05;
ampl_shift = 0.05;
ampl_num = 100;
ampls = linspace(ampl_start, ampl_start + (ampl_num - 1) * ampl_shift, ampl_num);
T_start = 0.5;
T_shift = 0.05;
T_num = 91;
Ts = linspace(T_start * (ps_drv_part_1 + ps_drv_part_2) * ps_prm_alpha, (T_start + (T_num - 1) * T_shift) * (ps_drv_part_1 + ps_drv_part_2) * ps_prm_alpha, T_num);
d = 1.0;
g = 10.0;

num_runs = 1;

alphas = zeros(ampl_num, T_num);
decades = zeros(ampl_num, T_num);

for T_id = 1:T_num
    
    T = T_start + T_shift * (T_id - 1);
    fprintf('T = %0.16e\n', T);
    
    for ampl_id = 1:ampl_num
        
        ampl = ampls(ampl_id);
        fprintf('ampl = %0.16e\n', ampl);
        
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
            
            path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/obs_%d_%d_%d/N_%d_%d/diss_%d_%0.4f/drv_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
                data_path, ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                ex_deep, ...
                rk_ns, ...
                num_tp_periods, ...
                num_obs_periods, ...
                num_random_obs, ...
                random_obs_seed, ...
                random_obs_type, ...
                ps_num_spins, ...
                ps_num_photons_states, ...
                diss_type, ...
                ps_diss_w, ...
                ps_drv_part_1 * T, ...
                ps_drv_part_2 * T, ...
                ampl, ...
                ps_prm_alpha, ...
                d, ...
                g, ...
                start_type, ...
                start_state, ...
                ss);
            
            suffix = sprintf("setup(%d_%d_%d)_rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                ss, ...
                mns, ...
                ps_num_spins, ...
                ps_num_photons_states, ...
                diss_type, ...
                ps_diss_w, ...
                ps_drv_part_1 * T, ...
                ps_drv_part_2 * T, ...
                ampl, ...
                ps_prm_alpha, ...
                d, ...
                g, ...
                start_type, ...
                start_state);
            
            for trajectory_id = 1:num_trajectories
                
                path = sprintf('%s/jump_times_%d_%s.txt', path_to_folder, trajectory_id - 1, suffix);
                data = importdata(path);
                path = sprintf('%s/diss_types_%d_%s.txt', path_to_folder, trajectory_id - 1, suffix);
                diss_types = importdata(path);
				
				if size(data, 1) > 1
				
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
            
        end
        
        non_inc_count = non_inc_count;
        norm = sum(curr_pdf);
        
        for bin_id = 1 : num_bins
            curr_pdf(bin_id) = curr_pdf(bin_id) / (norm * bin_diff(bin_id));
        end
        
        norm_check = 0;
        for bin_id = 1 : num_bins
            norm_check = norm_check + curr_pdf(bin_id) * bin_diff(bin_id);
        end
        norm_check = norm_check;
        
        [x_min, x_max] = oqs_find_range(bin_centers, curr_pdf);
        [alpha, coef, R2, yy, R2_, cnt] = oqs_powerlaw_regression(bin_centers, curr_pdf, x_min, x_max);
        left = find(bin_centers==x_min);
        right = find(bin_centers==x_max);
        
        if ~isnan(R2)
            alphas(ampl_id, T_id) = abs(alpha);
            decades(ampl_id, T_id) = log10(x_max) - log10(x_min);
        else
            alphas(ampl_id, T_id) = 0;
            decades(ampl_id, T_id) = 0;
        end
        
    end
    
end


suffix_save = sprintf("diss(%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_var)_prm(%0.4f_%0.4f_%0.4f)", ...
    diss, ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    diss_type, ...
    ps_diss_w, ...
    ps_drv_part_1, ...
    ps_drv_part_2, ...
    ps_prm_alpha, ...
    d, ...
    g);

fig = figure;
hLine = imagesc(ampls, Ts, alphas');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$T$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\alpha');
set(gca,'YDir','normal');
savefig(sprintf('%s/tbj_alphas_from_ampl_and_T_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/tbj_alphas_from_ampl_and_T_%s.pdf', home_figures_path, suffix_save));

close(fig);

fig = figure;

hLine = imagesc(ampls, Ts, decades');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$T$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'decades');
set(gca,'YDir','normal');
savefig(sprintf('%s/tbj_decades_from_ampl_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/tbj_decades_from_ampl_%s.pdf', home_figures_path, suffix_save));

close(fig);

