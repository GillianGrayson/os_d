clear all;

eps = 1.0e-8;

data_path = '../../../../data/cluster/unn/qjx';

sys_id = 0;
task_id = 1;
prop_id = 1;
seed = 0;
mns = 1000000;
num_trajectories = 32;
num_tp_periods = 100;
num_obs_periods = 1000;
ex_deep = 16;
rk_ns = 10000;

N = 100;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;
dimer_drv_type = 1; 
dimer_drv_ampl = 3.4; 
dimer_drv_freq = 1; 
dimer_drv_phase = 0; 
dimer_prm_E = 0;
dimer_prm_U = 0.01;
dimer_prm_J = 1; 
start_type = 0;
start_state = 0;

ss = 0;

sys_size = N + 1;

tr_id = 1;

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

suffix = sprintf('rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
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

fn = sprintf('%s/mean_evo_%s.txt', path_to_folder, suffix);
mean_evo = importdata(fn);
xs = mean_evo(:, tr_id);

fn = sprintf('%s/energy_evo_%s.txt', path_to_folder, suffix);
energy_evo = importdata(fn);
ys = energy_evo(:, tr_id);

num_bins_x = 400;
num_bins_y = 400;
x_bins = zeros(num_bins_x, 1);
y_bins = zeros(num_bins_y, 1);
z_data = zeros(num_bins_x, num_bins_y);

x_min = min(xs);
x_max = max(xs);
x_shift = (x_max - x_min) / num_bins_x;
x_bins = linspace(x_min + 0.5 * x_shift, x_max - 0.5 * x_shift, num_bins_x);
x_pdf = zeros(num_bins_x, 1);

y_min = min(ys);
y_max = max(ys);
y_shift = (y_max - y_min) / num_bins_y;
y_bins = linspace(y_min + 0.5 * y_shift, y_max - 0.5 * y_shift, num_bins_y);

for p_id = 1 : size(xs, 1)
    x = xs(p_id);
    y = ys(p_id);
    
    x_id = floor((x - x_min) / (x_max - x_min + eps) * num_bins_x) + 1;
    y_id = floor((y - y_min) / (y_max - y_min + eps) * num_bins_y) + 1;
    
    x_pdf(x_id) = x_pdf(x_id) + 1;
    
    z_data(x_id, y_id) = z_data(x_id, y_id) + 1;  
end

x_pdf = x_pdf / (size(xs, 1) * x_shift);
norm = sum(x_pdf) * x_shift
x_pdf = x_pdf / max(x_pdf);


z_data = z_data / (size(xs, 1) * x_shift * y_shift);
norm = sum(sum(z_data)) *  x_shift * y_shift


fig = figure;
hLine = imagesc(x_bins, y_bins, z_data');
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$e$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
% set(gca, 'XAxisLocation', 'top');
% set(gca,'xticklabel',[]);
% set(gca, 'Position', [0.15 0.50 0.70 0.40])
hold all;

propertyeditor(fig)


fig = figure;
hLine = plot(x_bins, x_pdf);
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
xlim([0 N])
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex');
hold all;
propertyeditor(fig)