clear all;

ss = 1;
qj_mns = 1000000;
N = 100;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;
prm_E = 0.0;
prm_U = 0.01;
prm_J = 1.0;
start_type = 0;
start_state = 0;
eps = 1.0e-8;

tr_id = 2;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
    ss, ...
    qj_mns, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase, ...
    prm_E, ...
    prm_U, ...
    prm_J, ...
    start_type, ...
    start_state);

fn = sprintf('%s/mean_evo_%s.txt', data_path, suffix);
mean_evo = importdata(fn);
xs = mean_evo(:, tr_id);

fn = sprintf('%s/energy_evo_%s.txt', data_path, suffix);
energy_evo = importdata(fn);
ys = energy_evo(:, tr_id);

num_bins_x = 100;
num_bins_y = 101;
x_bins = zeros(num_bins_x, 1);
y_bins = zeros(num_bins_y, 1);
z_data = zeros(num_bins_x, num_bins_y);

x_min = min(xs);
x_max = max(xs);
x_shift = (x_max - x_min) / num_bins_x;
x_bins = linspace(x_min + 0.5 * x_shift, x_max - 0.5 * x_shift, num_bins_x);

y_min = min(ys);
y_max = max(ys);
y_shift = (y_max - y_min) / num_bins_y;
y_bins = linspace(y_min + 0.5 * y_shift, y_max - 0.5 * y_shift, num_bins_y);

for p_id = 1 : size(xs, 1)
    x = xs(p_id);
    y = ys(p_id);
    
    x_id = floor((x - x_min) / (x_max - x_min + eps) * num_bins_x) + 1;
    y_id = floor((y - y_min) / (y_max - y_min + eps) * num_bins_y) + 1;
    
    z_data(x_id, y_id) = z_data(x_id, y_id) + 1;  
end

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