clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 0; 
task_id = 1;
prop_id = 1;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 100;
num_obs_periods = 100;

N = 100;

diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;

dimer_drv_type = 1;
dimer_drv_ampl = 3.4;
dimer_drv_freq = 1.0;
dimer_drv_phase = 0.0;

dimer_prm_E = 0.0;
dimer_prm_U = 0.15;
dimer_prm_J = 1.0;

start_type = 0;
start_state = 0;

ss = 0;
mns = 1000000;

num_trajectories = 10;
num_runs = 10;

sys_size = N + 1;

adaptive_axes = 0;
num_bins_x = sys_size;
num_bins_y = sys_size;

mean_evo = zeros(num_obs_periods + 1, num_trajectories * num_runs);
energy_evo = zeros(num_obs_periods + 1, num_trajectories * num_runs);

for run_id = 1:num_runs
    
    ss = (run_id - 1) * num_trajectories;
    
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
    mean_evo_curr = importdata(fn);
    mean_evo(:, ss + 1: ss + num_trajectories) = mean_evo_curr;
    
    fn = sprintf('%s/energy_evo_%s.txt', path_to_folder, suffix);
    energy_evo_curr = importdata(fn);
    energy_evo(:, ss + 1: ss + num_trajectories) = energy_evo_curr;
end

z_data = zeros(num_bins_x, num_bins_y);

x_min = min(min(mean_evo)) - 1e-4;
x_max = max(max(mean_evo)) + 1e-4;
x_shift = (x_max - x_min) / num_bins_x;
x_bins = linspace(x_min + 0.5 * x_shift, x_max - 0.5 * x_shift, num_bins_x);


y_min = min(min(energy_evo)) - 1e-4;
y_max = max(max(energy_evo)) + 1e-4;
y_shift = (y_max - y_min) / num_bins_y;
y_bins = linspace(y_min + 0.5 * y_shift, y_max - 0.5 * y_shift, num_bins_y);

total_num_trajectories = num_trajectories * num_runs;

for tr_id = 1:total_num_trajectories
    
    xs = mean_evo(:, tr_id);
    ys = energy_evo(:, tr_id);
    
    for p_id = 1 : size(xs, 1)
        x = xs(p_id);
        y = ys(p_id);

        x_id = floor((x - x_min) / (x_max - x_min + eps) * num_bins_x) + 1;
        y_id = floor((y - y_min) / (y_max - y_min + eps) * num_bins_y) + 1;

        z_data(x_id, y_id) = z_data(x_id, y_id) + 1;  
    end

end

z_data = z_data / (size(mean_evo, 1) * size(mean_evo, 2) * x_shift * y_shift);
norm = sum(sum(z_data)) *  x_shift * y_shift
z_data = z_data / max(max(z_data));

fig = figure;
hLine = imagesc(x_bins, y_bins, z_data');
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex')
set(gca, 'FontSize', 30);
ylabel('$e$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
hold all;

suffix_save = sprintf('N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
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

savefig(sprintf('%s/mean_energy_pdf_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/mean_energy_pdf_%s.pdf', home_figures_path, suffix_save));

