clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 3;
task_id = 7;
prop_id = 0;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 1000;
num_obs_periods = 1000;

Nc = 8;

W_seed_s = 1;
W_seed_num = 100;
W_seed_f = W_seed_s + W_seed_num - 1;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

U = 1;
J = 1;

start_type = 0;
start_state = 49;

lpn_delta_f = 0.001;
lpn_delta_s = 0.000001;

ss = 0;
mns = 1000000;

num_trajectories = 2000;
num_target_trajectories = num_trajectories / 2;
num_runs = 1;

lambda_min = -0.15;
lambda_max = 0.15;
lambda_num_bins = 200;
lambda_shift = (lambda_max - lambda_min) / lambda_num_bins;
lambda_bins = linspace(lambda_min + 0.5 * lambda_shift, lambda_max - 0.5 * lambda_shift, lambda_num_bins);

W_begin = 0.2;
W_step = 0.2;
W_num = 100;
Ws = zeros(W_num, 1);

lambdas_2d_pdf = zeros(W_num, lambda_num_bins);
lambdas_2d_pdf_deep = zeros(W_num, lambda_num_bins);

for W_id = 1:W_num
    
    W = W_begin + W_step * (W_id - 1)
    Ws(W_id) = W;
    
    lambdas_pdf = zeros(lambda_num_bins, 1);
    lambdas_pdf_deep = zeros(lambda_num_bins, 1);
    
    non_inc = 0;
    non_inc_deep = 0;
    
    for W_seed = W_seed_s : W_seed_f
        
        for run_id = 1:num_runs
            
            ss = (run_id - 1) * num_trajectories;
            
            path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/Nс_%d/rnd_%d_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/lpn_%0.4f_%0.4f/ss_%d', ...
                data_path, ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                ex_deep, ...
                rk_ns, ...
                num_tp_periods, ...
                num_obs_periods, ...
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
                log10(lpn_delta_s), ...
                log10(lpn_delta_f), ...
                ss);
            
            suffix = sprintf('rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%0.4f_%0.4f)', ...
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
                log10(lpn_delta_s), ...
                log10(lpn_delta_f));
            
            fn = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
            lambda_data = importdata(fn);
            curr_lambdas = lambda_data(num_target_trajectories + 1 : end);
            curr_lambdas_mean = mean(curr_lambdas);
            
            if curr_lambdas_mean < lambda_max && curr_lambdas_mean >= lambda_min
                curr_id = floor((curr_lambdas_mean - lambda_min) / (lambda_max - lambda_min + eps) * lambda_num_bins) + 1;
                lambdas_pdf(curr_id) = lambdas_pdf(curr_id) + 1;
            else
                non_inc = non_inc + 1;
            end
            
            for l_id = 1: size(curr_lambdas, 1)
                curr_l = curr_lambdas(l_id);
                if curr_l < lambda_max && curr_l >= lambda_min
                    curr_id = floor((curr_l - lambda_min) / (lambda_max - lambda_min + eps) * lambda_num_bins) + 1;
                    lambdas_pdf_deep(curr_id) = lambdas_pdf_deep(curr_id) + 1;
                else
                    non_inc_deep = non_inc_deep + 1;
                end
            end
            
            
        end 
    end
    
    non_inc = non_inc
    sum_lambda_pdf = sum(lambdas_pdf);
    lambdas_pdf = lambdas_pdf / (sum_lambda_pdf * lambda_shift);
    norm = sum(lambdas_pdf) * lambda_shift;
    norm_diff = 1.0 - norm
    lambdas_2d_pdf(W_id, :) = lambdas_pdf;
    
    non_inc_deep = non_inc_deep
    sum_lambda_pdf_deep = sum(lambdas_pdf_deep);
    lambdas_pdf_deep = lambdas_pdf_deep / (sum_lambda_pdf_deep * lambda_shift);
    norm = sum(lambdas_pdf_deep) * lambda_shift;
    norm_diff_deep = 1.0 - norm
    lambdas_2d_pdf_deep(W_id, :) = lambdas_pdf_deep;
    
end

for W_id = 1:W_num
    curr_max = max(lambdas_2d_pdf(W_id, :));
    lambdas_2d_pdf(W_id, :) = lambdas_2d_pdf(W_id, :) / curr_max;
    
    curr_max = max(lambdas_2d_pdf_deep(W_id, :));
    lambdas_2d_pdf_deep(W_id, :) = lambdas_2d_pdf_deep(W_id, :) / curr_max;
end

suffix_save =  sprintf('Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(var_%0.4f_%0.4f)_start(%d_%d)_lpn(%0.4f_%0.4f)', ...
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
    log10(lpn_delta_s), ...
    log10(lpn_delta_f));


fig = figure;
hLine = imagesc(Ws, lambda_bins, lambdas_2d_pdf');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '');
set(gca,'YDir','normal');

savefig(sprintf('%s/lambda_pdf_from_W_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_pdf_from_W_%s.pdf', home_figures_path, suffix_save));

close(fig)

fig = figure;
hLine = imagesc(Ws, lambda_bins, lambdas_2d_pdf_deep');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '');
set(gca,'YDir','normal');

savefig(sprintf('%s/lambda_pdf_deep_from_W_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_pdf_deep_from_W_%s.pdf', home_figures_path, suffix_save));

close(fig)
