clear all;

eps = 1.0e-8;

adaptive_axes = 0;

data_path = '../../../../data/cluster/unn/qjx';

sys_id = 0;
task_id = 1;
prop_id = 1;
seed = 0;
mns = 1000000;
num_trajectories = 32;
num_tp_periods = 1000;
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
dimer_prm_U = 0.1125;
dimer_prm_J = 1; 
start_type = 0;
start_state = 0;

sys_size = N + 1;

tr_id = 1;

num_runs = 20;
mean_evo = zeros(num_obs_periods + 1, num_trajectories * num_runs);
energy_evo = zeros(num_obs_periods + 1, num_trajectories * num_runs);

for run_id = 1:num_runs
    
    ss = 0 + (run_id - 1) * num_trajectories;

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
    % normalization
    mean_evo_curr = mean_evo_curr / N;
    mean_evo(:, ss + 1: ss + num_trajectories) = mean_evo_curr;
    
    fn = sprintf('%s/energy_evo_%s.txt', path_to_folder, suffix);
    energy_evo_curr = importdata(fn);
    energy_evo(:, ss + 1: ss + num_trajectories) = energy_evo_curr;

end

num_bins_x = 400;
num_bins_y = 400;
z_data = zeros(num_bins_x, num_bins_y);

if adaptive_axes == 1
    x_min = min(min(mean_evo));
    x_max = max(max(mean_evo));
else
    x_min = 0;
    x_max = 1;
end
x_shift = (x_max - x_min) / num_bins_x;
x_bins = linspace(x_min + 0.5 * x_shift, x_max - 0.5 * x_shift, num_bins_x);
x_pdf = zeros(num_bins_x, 1);

if adaptive_axes == 1
    y_min = min(min(energy_evo));
    y_max = max(max(energy_evo));
else
    y_min = min(min(energy_evo));
    y_max = max(max(energy_evo));
end
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

        x_pdf(x_id) = x_pdf(x_id) + 1;

        z_data(x_id, y_id) = z_data(x_id, y_id) + 1;  
    end

end

x_pdf = x_pdf / (size(mean_evo, 1) * size(mean_evo, 2) * x_shift);
norm = sum(x_pdf) * x_shift
x_pdf = x_pdf / max(x_pdf);


z_data = z_data / (size(mean_evo, 1) * size(mean_evo, 2) * x_shift * y_shift);
norm = sum(sum(z_data)) *  x_shift * y_shift
z_data = z_data / max(max(z_data));

fig = figure;
subplot(2,1,2);
hLine = plot(x_bins, x_pdf);
set(gca, 'Position', [0.15 0.15 0.70 0.25])
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
xlim([0 1])
set(gca, 'FontSize', 30);
ylabel('$|\psi_{n,n}|^2$', 'Interpreter', 'latex');
hold all;

subplot(2,1,1);
hLine = imagesc(x_bins, y_bins, z_data');
set(gca, 'FontSize', 30);
set(gca,'xticklabel',[]);
set(gca, 'FontSize', 30);
ylabel('$e$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
xlim([0 1])
set(gca, 'Position', [0.15 0.41 0.70 0.45])
hold all;

propertyeditor(fig)

suffix = sprintf('N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
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

fn = sprintf('mean_energy_pdf_%s', suffix);

h=gcf;
savefig(h, sprintf('%s.fig', fn))
set(h,'PaperOrientation','landscape'); 
set(h,'PaperUnits','normalized'); 
set(h,'PaperPosition', [0 0 1 1]); 
print(gcf, '-dpdf', sprintf('%s.pdf', fn));
