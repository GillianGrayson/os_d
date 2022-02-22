clear all;

home_figures_path = '/home/ivanchen/yusipov/os_d/figures';
data_path = '/data/condmat/ivanchen/yusipov/os_d/qjx';

sys_id = 3;
task_id = 7;
prop_id = 0;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 500;
num_obs_periods = 500;

T = 10;

Nc = 12;

num_rnd_obs = 1;
rnd_obs_seed = 100;
rnd_obs_type = 2;

W_seed_s = 1;
W_seed_num = 100;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

U = 1;
J = 1;

start_type = 0;
start_state = 0;

lpn_type = -1;
lpn_delta_f_h = 1e-6;
lpn_delta_s = 1e-6;
lpn_delta_f_l = 1e-6;

ss = 0;
mns = 1000000;

num_trajectories = 200;
num_target_trajectories = num_trajectories / 2;
num_runs = 1;

W_begin = 0.2;
W_step = 0.2;
W_num = 100;
Ws = zeros(W_num, 1);

all_lambdas = zeros(W_num, W_seed_num * num_runs * num_target_trajectories);
all_renorms = zeros(W_num, W_seed_num * num_runs * num_target_trajectories);

for W_id = 1:W_num
    
    W = W_begin + W_step * (W_id - 1)
    Ws(W_id) = W;
       
    for W_seed_id = 1 : W_seed_num
        
        W_seed = W_seed_s + W_seed_id - 1;
        
        for run_id = 1:num_runs
            
            ss = (run_id - 1) * num_trajectories;
            
            path_to_folder = sprintf('%s/main_%d_%d_%d/T_%0.4f/lpn_%d_%0.4f_%0.4f_%0.4f/run_%d_%d_%d_%d/Nс_%d/obs_%d_%d_%d/rnd_%d_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
				data_path, ...
				sys_id, ...
				task_id, ...
				prop_id, ...
				T, ...
				lpn_type, ...
				log10(lpn_delta_s), ...
				log10(lpn_delta_f_h), ...
				log10(lpn_delta_f_l), ...
				ex_deep, ...
				rk_ns, ...
				num_tp_periods, ...
				num_obs_periods, ...
				Nc, ...
				num_rnd_obs, ...
				rnd_obs_seed, ...
				rnd_obs_type, ...
				W_seed, ...
				W_mns, ...
				diss_type, ...
				diss_phase, ...
				diss_gamma, ...
				W, ...
				U, ...
				J, ...
				start_type, ...
				start_state, ...
				ss);
            
            suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                ss, ...
                mns, ...
                Nc, ...
                W_seed, ...
                W_mns, ...
                diss_type, ...
                diss_phase, ...
                diss_gamma, ...
                W, ...
                U, ...
                J, ...
                start_type, ...
                start_state, ...
                lpn_type, ...
                log10(lpn_delta_s), ...
                log10(lpn_delta_f_h), ...
                log10(lpn_delta_f_l));
            s_index = (W_seed_id - 1) * num_runs * num_target_trajectories + (num_runs - 1) * num_target_trajectories + 1;
            f_index = (W_seed_id - 1) * num_runs * num_target_trajectories + (num_runs - 1) * num_target_trajectories + num_target_trajectories;
            
            fn = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
            curr_data = importdata(fn);
            curr_lambdas = curr_data(num_target_trajectories + 1 : end);
            all_lambdas(W_id, s_index : f_index) = curr_lambdas;
            
            fn = sprintf('%s/num_renorms_%s.txt', path_to_folder, suffix);
            curr_data = importdata(fn);
            curr_renorms = curr_data(num_target_trajectories + 1 : end);
            all_renorms(W_id, s_index : f_index) = curr_renorms;
        end 
    end
end

lambda_min = min(all_lambdas, [], 'all') - 1e-10;
lambda_max = max(all_lambdas, [], 'all') + 1e-10;
lambda_num_bins = 400;
lambda_shift = (lambda_max - lambda_min) / lambda_num_bins;
lambda_bins = linspace(lambda_min + 0.5 * lambda_shift, lambda_max - 0.5 * lambda_shift, lambda_num_bins);

lambdas_2d_pdf_deep = zeros(W_num, lambda_num_bins);
lambdas_2d_pdf = zeros(W_num, lambda_num_bins);
lambdas = zeros(W_num, 1);

renorm_min = min(all_renorms, [], 'all') - 1e-10;
renorm_max = max(all_renorms, [], 'all') + 1e-10;
renorm_num_bins = 400;
renorm_shift = (renorm_max - renorm_min) / renorm_num_bins;
renorm_bins = linspace(renorm_min + 0.5 * renorm_shift, renorm_max - 0.5 * renorm_shift, renorm_num_bins);

renorms_2d_pdf_deep = zeros(W_num, renorm_num_bins);
renorms_2d_pdf = zeros(W_num, renorm_num_bins);
renorms = zeros(W_num, 1);

for W_id = 1:W_num
    
    W = W_begin + W_step * (W_id - 1);
    
    lambdas_pdf_deep = zeros(lambda_num_bins, 1);
    lambdas_pdf = zeros(lambda_num_bins, 1);
    avg_lambdas = 0;
    
    renorms_pdf_deep = zeros(renorm_num_bins, 1);
    renorms_pdf = zeros(renorm_num_bins, 1);
    avg_renorms = 0;
    
    for W_seed_id = 1 : W_seed_num
        
        W_seed = W_seed_s + W_seed_id - 1;
        
        for run_id = 1:num_runs
            
            ss = (run_id - 1) * num_trajectories;
            
            s_index = (W_seed_id - 1) * num_runs * num_trajectories + (num_runs - 1) * num_trajectories + 1;
            f_index = (W_seed_id - 1) * num_runs * num_trajectories + (num_runs - 1) * num_trajectories + num_target_trajectories;
            
            curr_lambdas = all_lambdas(W_id, s_index : f_index);
            for l_id = 1: size(curr_lambdas, 1)
                curr_l = curr_lambdas(l_id);
                curr_id = floor((curr_l - lambda_min) / (lambda_max - lambda_min) * lambda_num_bins) + 1;
                lambdas_pdf_deep(curr_id) = lambdas_pdf_deep(curr_id) + 1;
            end
            
            curr_lambdas_mean = mean(curr_lambdas);
            curr_id = floor((curr_lambdas_mean - lambda_min) / (lambda_max - lambda_min ) * lambda_num_bins) + 1;
            lambdas_pdf(curr_id) = lambdas_pdf(curr_id) + 1;

			avg_lambdas = avg_lambdas + curr_lambdas_mean;
            
            
            curr_renorms = all_renorms(W_id, s_index : f_index);
            for l_id = 1: size(curr_renorms, 1)
                curr_l = curr_renorms(l_id);
                curr_id = floor((curr_l - renorm_min) / (renorm_max - renorm_min) * renorm_num_bins) + 1;
                renorms_pdf_deep(curr_id) = renorms_pdf_deep(curr_id) + 1;
            end
            
            curr_renorms_mean = mean(curr_renorms);
            curr_id = floor((curr_renorms_mean - renorm_min) / (renorm_max - renorm_min ) * renorm_num_bins) + 1;
            renorms_pdf(curr_id) = renorms_pdf(curr_id) + 1;

			avg_renorms = avg_renorms + curr_renorms_mean;
        end 
    end
    
    sum_lambda_pdf_deep = sum(lambdas_pdf_deep);
    lambdas_pdf_deep = lambdas_pdf_deep / (sum_lambda_pdf_deep * lambda_shift);
    norm = sum(lambdas_pdf_deep) * lambda_shift;
    norm_diff_deep = 1.0 - norm
    lambdas_2d_pdf_deep(W_id, :) = lambdas_pdf_deep;
    curr_max = max(lambdas_2d_pdf_deep(W_id, :));
    lambdas_2d_pdf_deep(W_id, :) = lambdas_2d_pdf_deep(W_id, :) / curr_max;
    
    sum_lambda_pdf = sum(lambdas_pdf);
    lambdas_pdf = lambdas_pdf / (sum_lambda_pdf * lambda_shift);
    norm = sum(lambdas_pdf) * lambda_shift;
    norm_diff = 1.0 - norm
    lambdas_2d_pdf(W_id, :) = lambdas_pdf;  
    curr_max = max(lambdas_2d_pdf(W_id, :));
    lambdas_2d_pdf(W_id, :) = lambdas_2d_pdf(W_id, :) / curr_max;
	
	avg_lambdas = avg_lambdas / (W_seed_num * num_runs)
	lambdas(W_id) = avg_lambdas; 
    
    
    sum_renorm_pdf_deep = sum(renorms_pdf_deep);
    renorms_pdf_deep = renorms_pdf_deep / (sum_renorm_pdf_deep * renorm_shift);
    norm = sum(renorms_pdf_deep) * renorm_shift;
    norm_diff_deep = 1.0 - norm
    renorms_2d_pdf_deep(W_id, :) = renorms_pdf_deep;
    curr_max = max(renorms_2d_pdf_deep(W_id, :));
    renorms_2d_pdf_deep(W_id, :) = renorms_2d_pdf_deep(W_id, :) / curr_max;
    
    sum_renorm_pdf = sum(renorms_pdf);
    renorms_pdf = renorms_pdf / (sum_renorm_pdf * renorm_shift);
    norm = sum(renorms_pdf) * renorm_shift;
    norm_diff = 1.0 - norm
    renorms_2d_pdf(W_id, :) = renorms_pdf;  
    curr_max = max(renorms_2d_pdf(W_id, :));
    renorms_2d_pdf(W_id, :) = renorms_2d_pdf(W_id, :) / curr_max;
	
	avg_renorms = avg_renorms / (W_seed_num * num_runs)
	renorms(W_id) = avg_renorms;   
end

suffix_save =  sprintf('T(%0.4f)_setup(%d_%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(var_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    T, ...
	sys_id, ...
    task_id, ...
    prop_id, ...
    Nc, ...
    W_seed, ...
    W_mns, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    U, ...
    J, ...
    start_type, ...
    start_state, ...
	lpn_type, ...
    log10(lpn_delta_s), ...
    log10(lpn_delta_f_h), ...
    log10(lpn_delta_f_l));

fig = figure;
hLine = imagesc(Ws, lambda_bins, lambdas_2d_pdf_deep');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\lambda_d$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF', 'Interpreter', 'latex');
set(gca,'YDir','normal');
hold all;
hLine = plot(Ws, lambdas);
savefig(sprintf('%s/lambda_deep_pdf_from_W_%s.fig', home_figures_path, suffix_save));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_deep_pdf_from_W_%s.pdf', home_figures_path, suffix_save));
close(fig)

fig = figure;
hLine = imagesc(Ws, lambda_bins, lambdas_2d_pdf');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF', 'Interpreter', 'latex');
set(gca,'YDir','normal');
hold all;
hLine = plot(Ws, lambdas);
savefig(sprintf('%s/lambda_pdf_from_W_rnd_obs_%s.fig', home_figures_path, suffix_save));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_pdf_from_W_rnd_obs_%s.pdf', home_figures_path, suffix_save));
close(fig)

fig = figure;
hLine = imagesc(Ws, renorm_bins, renorms_2d_pdf_deep');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$renorms$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF', 'Interpreter', 'latex');
set(gca,'YDir','normal');
hold all;
hLine = plot(Ws, renorms);
savefig(sprintf('%s/renorm_deep_pdf_from_W_%s.fig', home_figures_path, suffix_save));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/renorm_deep_pdf_from_W_%s.pdf', home_figures_path, suffix_save));
close(fig)

fig = figure;
hLine = imagesc(Ws, renorm_bins, renorms_2d_pdf');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$renorms$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF', 'Interpreter', 'latex');
set(gca,'YDir','normal');
hold all;
hLine = plot(Ws, renorms);
savefig(sprintf('%s/renorm_pdf_from_W_rnd_obs_%s.fig', home_figures_path, suffix_save));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/renorm_pdf_from_W_rnd_obs_%s.pdf', home_figures_path, suffix_save));
close(fig)