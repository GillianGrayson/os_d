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
dimer_prm_U = 0.25;
dimer_prm_J = 1.0;

start_type = 0;
start_state = 0;

ss = 0;
mns = 1000000;

num_trajectories = 10;
num_runs = 10;

sys_size = N + 1;

num_bins_x = 100;

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

mean_x = mean(mean_evo(:));
mean_y = mean(energy_evo(:));

total_num_trajectories = num_trajectories * num_runs;

phase_evo = zeros(num_obs_periods + 1, total_num_trajectories);

for tr_id = 1:total_num_trajectories
    
    xs = mean_evo(:, tr_id);
    ys = energy_evo(:, tr_id);
    
    for p_id = 1 : size(xs, 1)
        x = xs(p_id) - mean_x;
        y = ys(p_id) - mean_y;
        
        if x > 0 && y >= 0
            phase_evo(p_id, tr_id) = atan(y/x);
        elseif x > 0 && y < 0
            phase_evo(p_id, tr_id) = atan(y/x) + 2 * pi;
        elseif x < 0
            phase_evo(p_id, tr_id) = atan(y/x) + pi;
        end
    end
end

phase_difference = zeros(num_obs_periods, total_num_trajectories);
for tr_id = 1:total_num_trajectories
    phase_difference(:, tr_id) = diff(phase_evo(:, tr_id)) / (2 * pi);
end

x_min = min(phase_difference(:)) - 1e-3;
x_max = max(phase_difference(:)) + 1e-3;
x_shift = (x_max - x_min) / num_bins_x;
x_bins = linspace(x_min + 0.5 * x_shift, x_max - 0.5 * x_shift, num_bins_x);

x_pdf = zeros(num_bins_x, 1);

for tr_id = 1:total_num_trajectories
    
    xs = phase_difference(:, tr_id);
    
    for p_id = 1 : size(xs, 1)
        x = xs(p_id);
        x_id = floor((x - x_min) / (x_max - x_min + eps) * num_bins_x) + 1;
        x_pdf(x_id) = x_pdf(x_id) + 1;
    end
end
  
x_pdf = x_pdf / (num_obs_periods * total_num_trajectories * x_shift);

norm = sum(x_pdf) * x_shift;
norm_diff = 1.0 - norm

fig = figure;
hLine = plot(x_bins, x_pdf);
set(gca, 'FontSize', 30);
xlabel('phase diff', 'Interpreter', 'latex');
xlim([x_bins(1) x_bins(end)])
set(gca, 'FontSize', 30);
ylabel('$PDF$', 'Interpreter', 'latex');

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

savefig(sprintf('%s/phase_diff_pdf_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/phase_diff_pdf_%s.pdf', home_figures_path, suffix_save));

