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

W = 2;
U = 1;
J = 1;

start_type = 0;
start_state = 49;

lpn_type = 0;
lpn_delta_f_h = 1e-4;
lpn_delta_f_l = 1e-8;
lpn_delta_s = 1e-6;

ss = 0;
mns = 1000000;

num_trajectories = 200;
num_target_trajectories = num_trajectories / 2;
num_runs = 1;

obs_y_num_bins = 400;

non_inc_deep = 0;
all_obs_x = zeros(num_target_trajectories, W_seed_num * num_runs);
all_obs_y = zeros(num_target_trajectories, W_seed_num * num_runs);

W_seed_id = 1;
for W_seed = W_seed_s : W_seed_f
    
    for run_id = 1:num_runs
        
        curr_obs_id = (W_seed_id - 1) * num_runs + run_id;
        
        ss = (run_id - 1) * num_trajectories;
        
        path_to_folder = sprintf('%s/main_%d_%d_%d/lpn_%d_%0.4f_%0.4f_%0.4f/run_%d_%d_%d_%d/Nс_%d/rnd_%d_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
            data_path, ...
            sys_id, ...
            task_id, ...
            prop_id, ...
            lpn_type, ...
            log10(lpn_delta_s), ...
            log10(lpn_delta_f_h), ...
            log10(lpn_delta_f_l), ...
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
        
        fn = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
        obs_1_data = importdata(fn);
        curr_obs_1 = obs_1_data(num_target_trajectories + 1 : end);
        all_obs_y(:, curr_obs_id) = curr_obs_1;
        
        fn = sprintf('%s/num_renorms_%s.txt', path_to_folder, suffix);
        obs_2_data = importdata(fn);
        curr_obs_2 = obs_2_data(num_target_trajectories + 1 : end);
        all_obs_x(:, curr_obs_id) = curr_obs_2;
    end
    
    W_seed_id = W_seed_id + 1;
end

x_min = round(min(min(all_obs_x))) - 0.5;
x_max = round(max(max(all_obs_x))) + 0.5;
obs_x_num_bins = round(x_max - x_min);
x_shift = (x_max - x_min) / obs_x_num_bins;
x_bins = linspace(x_min + 0.5 * x_shift, x_max - 0.5 * x_shift, obs_x_num_bins);

y_min = min(min(all_obs_y)) - 1e-8;
y_max = max(max(all_obs_y)) + 1e-8;
y_shift = (y_max - y_min) / obs_y_num_bins;
y_bins = linspace(y_min + 0.5 * y_shift, y_max - 0.5 * y_shift, obs_y_num_bins);

pdf_deep_2d = zeros(obs_x_num_bins, obs_y_num_bins);
pdf_deep_x = zeros(obs_x_num_bins, 1);
pdf_deep_y = zeros(obs_y_num_bins, 1);
for i = 1 : num_target_trajectories
    for j = 1 : W_seed_num * num_runs
        curr_x = all_obs_x(i, j);
        curr_y = all_obs_y(i, j);
        
        x_id = floor((curr_x - x_min) / (x_max - x_min) * obs_x_num_bins) + 1;
        y_id = floor((curr_y - y_min) / (y_max - y_min) * obs_y_num_bins) + 1;
        
        pdf_deep_2d(x_id, y_id) = pdf_deep_2d(x_id, y_id) + 1;
        pdf_deep_x(x_id) = pdf_deep_x(x_id) + 1;
        pdf_deep_y(y_id) = pdf_deep_y(y_id) + 1;
    end
end

sum_pdf_deep_2d = sum(sum(pdf_deep_2d));
pdf_deep_2d = pdf_deep_2d / (sum_pdf_deep_2d * x_shift * y_shift);
norm = sum(sum(pdf_deep_2d)) * x_shift * y_shift;
norm_diff_deep_2d = 1.0 - norm

sum_pdf_deep_x = sum(pdf_deep_x);
pdf_deep_x = pdf_deep_x / (sum_pdf_deep_x * x_shift);
norm = sum(pdf_deep_x) * x_shift;
norm_diff_deep_x = 1.0 - norm

sum_pdf_deep_y = sum(pdf_deep_y);
pdf_deep_y = pdf_deep_y / (sum_pdf_deep_y * y_shift);
norm = sum(pdf_deep_y) * y_shift;
norm_diff_deep_y = 1.0 - norm

suffix_save = sprintf('setup(%d_%d_%d)_Nc(%d)_W_seed(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    Nc, ...
    W_seed_s, ...
    W_seed_f, ...
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

fig = figure;
hLine = imagesc(x_bins, y_bins, pdf_deep_2d');
set(gca, 'FontSize', 30);
xlabel('num renorms', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF', 'Interpreter', 'latex');
set(gca,'YDir','normal');

savefig(sprintf('%s/num_renorms_deep_lambda_deep_pdf_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/num_renorms_deep_lambda_deep_pdf_%s.pdf', home_figures_path, suffix_save));

close(fig)


fig = figure;
hLine = plot(x_bins, pdf_deep_x, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('num renorms', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex')

savefig(sprintf('%s/num_renorms_deep_pdf_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/num_renorms_deep_pdf_%s.pdf', home_figures_path, suffix_save));

close(fig)


fig = figure;
hLine = plot(y_bins, pdf_deep_y, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$\lambda$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex')

savefig(sprintf('%s/lambda_pdf_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_pdf_%s.pdf', home_figures_path, suffix_save));

close(fig)