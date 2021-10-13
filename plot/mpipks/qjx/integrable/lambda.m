clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 5;
task_id = 7;
prop_id = 0;


lpn_type = 0;
lpn_delta = 1e-6;
lpn_delta_f_h = lpn_delta;
lpn_delta_s = lpn_delta;
lpn_delta_f_l = lpn_delta;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 100;
num_obs_periods = 100;

N = 7; 
tau = 1;
k = -1;
T = 1.0;

num_seeds = 50;

ss = 0;
mns = 1000000;

num_trajectories = 100;
num_target_trajectories = num_trajectories / 2;

lambdas_glob = zeros(num_seeds, num_target_trajectories);

for seed = 1:num_seeds

	seed = seed

    path_to_folder = sprintf('%s/main_%d_%d_%d/T_%0.4f/lpn_%d_%0.4f/run_%d_%d_%d_%d/N_%d_tau_%d_k_%d/seed_%d/ss_%d', ...
        data_path, ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        T, ...
        lpn_type, ...
        log10(lpn_delta), ...
        ex_deep, ...
        rk_ns, ...
        num_tp_periods, ...
        num_obs_periods, ...
        N, ...
        tau, ...
        k, ...
        seed, ...
        ss);

    suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_N(%d)_seed(%d)_tau(%d)_k(%d)_T(%0.4f)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        ss, ...
        mns, ...
        N, ...
        seed, ...
        tau, ...
        k, ...
        T, ...
        lpn_type, ...
        log10(lpn_delta_s), ...
        log10(lpn_delta_f_h), ...
        log10(lpn_delta_f_l));
    
    
    fn = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
    lambdas_data = importdata(fn);
    lambdas = lambdas_data(num_target_trajectories + 1 : end);
    lambdas_glob(seed, :) = lambdas;
    
end

lambdas_glob_1d = lambdas_glob(:);
lambdas_glob_1d_mean = mean(lambdas_glob, 1);
suffix_save = sprintf('N(%d)_numSeeds(%d)_tau(%d)_k(%d)_T(%0.4f)_lpn(%d_%0.4f)', ...
    N, ...
    num_seeds, ...
    tau, ...
    k, ...
    T, ...
    lpn_type, ...
    log10(lpn_delta));
x_pos = 1;
    

fig = figure;
b = boxplot(lambdas_glob_1d, 'Notch', 'off', 'positions', x_pos, 'Colors', 'r');
set(gca, 'FontSize', 40);
all_items = handle(b);
tags = get(all_items,'tag');
idx = strcmpi(tags,'box');
boxes = all_items(idx);
set(all_items,'linewidth',3)
idx = strcmpi(tags,'Outliers');
outliers = all_items(idx);
set(outliers, 'visible', 'off')
hold all;
%xs = x_pos * ones(size(lambdas_glob_1d, 1), 1) + ((rand(size(lambdas_glob_1d))-0.5)/10);
%h = scatter(xs, lambdas_glob_1d, 100, 'o', 'LineWidth',  1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerEdgeAlpha', 0.6, 'MarkerFaceAlpha', 0.6);
%h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold all;
savefig(sprintf('%s/lambda_%s.fig', home_figures_path, suffix_save));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_%s.pdf', home_figures_path, suffix_save));
close(fig)

fig = figure;
b = boxplot(lambdas_glob_1d_mean, 'Notch', 'off', 'positions', x_pos, 'Colors', 'r');
set(gca, 'FontSize', 40);
all_items = handle(b);
tags = get(all_items,'tag');
idx = strcmpi(tags,'box');
boxes = all_items(idx);
set(all_items,'linewidth',3)
idx = strcmpi(tags,'Outliers');
outliers = all_items(idx);
set(outliers, 'visible', 'off')
hold all;
%xs = x_pos * ones(size(lambdas_glob_1d, 1), 1) + ((rand(size(lambdas_glob_1d))-0.5)/10);
%h = scatter(xs, lambdas_glob_1d, 100, 'o', 'LineWidth',  1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'red', 'MarkerEdgeAlpha', 0.6, 'MarkerFaceAlpha', 0.6);
%h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold all;
savefig(sprintf('%s/lambda_mean_%s.fig', home_figures_path, suffix_save));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_mean_%s.pdf', home_figures_path, suffix_save));
close(fig)

xlswrite(sprintf('%s/lambda_%s.xlsx', home_figures_path, suffix_save), lambdas_glob)

