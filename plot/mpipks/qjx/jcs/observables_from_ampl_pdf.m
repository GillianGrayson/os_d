clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';

data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

obs_bin_num = 500;

num_runs = 1;
num_trajectories = 200;


T = 1;

ampl_begin = 0.05;
ampl_step = 0.05;
ampl_num = 100;
ampls = zeros(ampl_num, 1);
lambdas = zeros(ampl_num, 1);

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
num_dumps = num_obs_periods + 1;

num_random_obs = 1;
random_obs_seed = 100;
random_obs_type = 2;

N = 220;

diss_type = 0;

jcs_drv_part_1 = 1.0 * T;
jcs_drv_part_2 = 1.0 * T;

jcs_prm_alpha = 5;

start_type = 0;
start_state = 0;

seed = 0;
mns = 1000000;

ampl_begin = 0.05;
ampl_step = 0.05;
ampl_num = 100;
ampls = zeros(ampl_num, 1);

num_trajectories_obs = num_trajectories;
beg_tr_id = 1;
end_tr_id = num_trajectories;
if task_id == 7
    num_trajectories_obs = num_trajectories / 2;
    beg_tr_id = 1;
    end_tr_id = num_trajectories / 2;
end

obs_real_all = zeros(ampl_num, num_runs * num_trajectories_obs, num_dumps);
obs_imag_all = zeros(ampl_num, num_runs * num_trajectories_obs, num_dumps);
obs_phot_all = zeros(ampl_num, num_runs * num_trajectories_obs, num_dumps);

for ampl_id = 1:ampl_num
    
    ampl_id = ampl_id
    ampl = ampl_begin + ampl_step * (ampl_id - 1);
    ampls(ampl_id) = ampl;
    
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
            jcs_drv_part_1, ...
            jcs_drv_part_2, ...
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
            jcs_drv_part_1, ...
            jcs_drv_part_2, ...
            ampl, ...
            jcs_prm_alpha, ...
            start_type, ...
            start_state, ...
            lpn_type, ...
            lpn_delta_s, ...
            lpn_delta_f_high, ...
            lpn_delta_f_low);
        
        
        path = sprintf('%s/spec_evo_%s.txt', path_to_folder, suffix);
        data_obs = importdata(path);
        
        path = sprintf('%s/mean_evo_%s.txt', path_to_folder, suffix);
        data_phot = importdata(path);
        
        
        for tr_id = beg_tr_id : end_tr_id
            
            tr_id_global = (run_id - 1) * num_trajectories_obs + tr_id;
            
            for d_id = 1 : size(data_obs, 1)
                
                curr_real = data_obs(d_id, (tr_id - 1) * 2 + 1);
                curr_imag = data_obs(d_id, (tr_id - 1) * 2 + 2);
                obs_real_all(ampl_id, tr_id_global, d_id) = curr_real;
                obs_imag_all(ampl_id, tr_id_global, d_id) = curr_imag;
                
                curr_phot = data_phot(d_id, tr_id);
                obs_phot_all(ampl_id, tr_id_global, d_id) = curr_phot;
            end
        end
    end
end

obs_real_all = obs_real_all * jcs_prm_alpha;
obs_imag_all = obs_imag_all * jcs_prm_alpha;

max_obs_real = max(obs_real_all, [], 'all') + 1e-6
min_obs_real = min(obs_real_all, [], 'all') - 1e-6
max_obs_imag = max(obs_imag_all, [], 'all') + 1e-6
min_obs_imag = min(obs_imag_all, [], 'all') - 1e-6

max_obs_phot = 250
min_obs_phot = 0

obs_shift_real = (max_obs_real - min_obs_real) / obs_bin_num;
obs_shift_imag = (max_obs_imag - min_obs_imag) / obs_bin_num;
obs_shift_phot = (max_obs_phot - min_obs_phot) / obs_bin_num;
obs_bins_real = zeros(obs_bin_num, 1);
obs_bins_imag = zeros(obs_bin_num, 1);
obs_bins_phot = zeros(obs_bin_num, 1);
for bin_id = 1 : obs_bin_num
    obs_bins_real(bin_id) = min_obs_real + ((bin_id - 1) + 0.5) * obs_shift_real;
    obs_bins_imag(bin_id) = min_obs_imag + ((bin_id - 1) + 0.5) * obs_shift_imag;
    obs_bins_phot(bin_id) = min_obs_phot + ((bin_id - 1) + 0.5) * obs_shift_phot;
end


pdf_obs_real = zeros(ampl_num, obs_bin_num);
pdf_obs_imag = zeros(ampl_num, obs_bin_num);
pdf_obs_phot = zeros(ampl_num, obs_bin_num);

for ampl_id = 1 : ampl_num
    
    for tr_id = 1 : num_runs * num_trajectories_obs
        
        for dump_id = 1 : num_dumps
            
            curr_real = obs_real_all(ampl_id, tr_id, dump_id);
            curr_imag = obs_imag_all(ampl_id, tr_id, dump_id);
            curr_phot = obs_phot_all(ampl_id, tr_id, dump_id);
            
            real_bin_id = floor((curr_real - min_obs_real) / obs_shift_real) + 1;
            imag_bin_id = floor((curr_imag - min_obs_imag) / obs_shift_imag) + 1;
            phot_bin_id = floor((curr_phot - min_obs_phot) / obs_shift_phot) + 1;
			
            pdf_obs_real(ampl_id, real_bin_id) = pdf_obs_real(ampl_id, real_bin_id) + 1;
            pdf_obs_imag(ampl_id, imag_bin_id) = pdf_obs_imag(ampl_id, imag_bin_id) + 1;
            pdf_obs_phot(ampl_id, phot_bin_id) = pdf_obs_phot(ampl_id, phot_bin_id) + 1;
            
        end
    end
    
    sum_real = sum(pdf_obs_real(ampl_id, :));
    pdf_obs_real(ampl_id, :) = pdf_obs_real(ampl_id, :) / (sum_real * obs_shift_real);
    norm_real_check = 1.0 - sum(pdf_obs_real(ampl_id, :)) * obs_shift_real;
	if abs(norm_real_check) > 1e-15
		norm_real_check = norm_real_check
	end
    pdf_obs_real(ampl_id, :) = pdf_obs_real(ampl_id, :) / max(pdf_obs_real(ampl_id, :));
    
    
    sum_imag = sum(pdf_obs_imag(ampl_id, :));
    pdf_obs_imag(ampl_id, :) = pdf_obs_imag(ampl_id, :) / (sum_imag * obs_shift_imag);
    norm_imag_check = 1.0 - sum(pdf_obs_imag(ampl_id, :)) * obs_shift_imag;
	if abs(norm_imag_check) > 1e-15
		norm_imag_check = norm_imag_check
	end
    pdf_obs_imag(ampl_id, :) = pdf_obs_imag(ampl_id, :) / max(pdf_obs_imag(ampl_id, :));
    
        
    sum_phot = sum(pdf_obs_phot(ampl_id, :));
    pdf_obs_phot(ampl_id, :) = pdf_obs_phot(ampl_id, :) / (sum_phot * obs_shift_phot);
    norm_phot_check = 1.0 - sum(pdf_obs_phot(ampl_id, :)) * obs_shift_phot;
	if abs(norm_phot_check) > 1e-15
		norm_phot_check = norm_phot_check
	end
    pdf_obs_phot(ampl_id, :) = pdf_obs_phot(ampl_id, :) / max(pdf_obs_phot(ampl_id, :));
    
end

suffix = sprintf('lpn(%d_%0.4f_%0.4f_%0.4f)_N(%d)_drv(%0.4f_%0.4f_var)_prm(%0.4f)_start(%d_%d)_time(%d_%d_%d_%d)_numTraj(%d)', ...
    lpn_type, ...
    lpn_delta_s, ...
    lpn_delta_f_high, ...
    lpn_delta_f_low,...
    N, ...
    jcs_drv_part_1, ...
    jcs_drv_part_2, ...
    jcs_prm_alpha, ...
    start_type, ...
    start_state, ...
    num_tp_periods, ...
    num_obs_periods, ...
    num_trajectories_obs);

fig = figure;
hLine = imagesc(ampls, obs_bins_real, pdf_obs_real');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$Re(\xi)$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF');
set(gca,'YDir','normal');

savefig(sprintf('%s/real_spec_from_ampl_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/real_spec_from_ampl_%s.pdf', home_figures_path, suffix));

close(fig);

fig = figure;
hLine = imagesc(ampls, obs_bins_real, pdf_obs_imag');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$Im(\xi)$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF');
set(gca,'YDir','normal');

savefig(sprintf('%s/imag_spec_from_ampl_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/imag_spec_from_ampl_%s.pdf', home_figures_path, suffix));

close(fig);


fig = figure;
hLine = imagesc(ampls, obs_bins_phot, pdf_obs_phot');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF');
set(gca,'YDir','normal');

savefig(sprintf('%s/mean_from_ampl_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/mean_from_ampl_%s.pdf', home_figures_path, suffix));

close(fig);