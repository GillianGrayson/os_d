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

N = 250;

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

num_target_trajectories = 1;

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

min_x = min(mean_evo(:));
max_x = max(mean_evo(:));
mean_evo = (mean_evo - min_x) / (max_x - min_x) * 2 - 1;
mean_x = mean(mean_evo(:))
mean_evo = mean_evo - mean_x;
min_x = min(mean_evo(:))
max_x = max(mean_evo(:))
	
min_y = min(energy_evo(:));
max_y = max(energy_evo(:));
energy_evo = (energy_evo - min_y) / (max_y - min_y) * 2 - 1;
mean_y = mean(energy_evo(:))
energy_evo = energy_evo - mean_y;
min_y = min(energy_evo(:))
max_y = max(energy_evo(:))

phase_diff = zeros(num_obs_periods, num_target_trajectories);

fig = figure;
for tr_id = 1:num_target_trajectories
    
    xs = mean_evo(:, tr_id);
    ys = energy_evo(:, tr_id);
    
    for p_id = 1 : size(xs, 1) - 1
        x1 = xs(p_id);
        x2 = xs(p_id + 1);
        y1 = ys(p_id);
        y2 = ys(p_id + 1);
        
        if x1 > 0 && y1 >= 0
            phase_1 = atan(y1/x1);
        elseif x1 > 0 && y1 < 0
            phase_1 = atan(y1/x1) + 2 * pi;
        elseif x1 < 0
            phase_1 = atan(y1/x1) + pi;
        end
		phase_1 = phase_1 / (2 * pi);
        
        if x2 > 0 && y2 >= 0
            phase_2 = atan(y2/x2);
        elseif x2 > 0 && y2 < 0
            phase_2 = atan(y2/x2) + 2 * pi;
        elseif x2 < 0
            phase_2 = atan(y2/x2) + pi;
        end
		phase_2 = phase_2 / (2 * pi);
        
		phase_diff_tmp = phase_2 - phase_1;
		if phase_diff_tmp < 0
			phase_diff_tmp = phase_diff_tmp + 1;
		end
		
        phase_diff(p_id, tr_id) = phase_diff_tmp;
        
        h_line = plot([x1 x2], [y1 y2]);
		
        legend(h_line, sprintf('%0.4f %0.4f %0.4f', phase_diff_tmp, phase_2, phase_1));
        hold all;
    end
end

hold all;
plot([min_x, max_x], [0 0], 'k', 'LineWidth', 3);
xlim([min_x, max_x])
hold all;
plot([0, 0], [min_y, max_y], 'k', 'LineWidth', 3);
ylim([min_y, max_y])

set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex')
set(gca, 'FontSize', 30);
ylabel('$e$', 'Interpreter', 'latex');
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

savefig(sprintf('%s/mean_energy_trajectories_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/mean_energy_trajectories_%s.pdf', home_figures_path, suffix_save));

phase_diff_min = 0;
phase_diff_max = 1;
phase_diff_num = 100;
phase_diff_shift = (phase_diff_max - phase_diff_min) / phase_diff_num;
phase_diff_bins = linspace(phase_diff_min + 0.5 * phase_diff_shift, phase_diff_max - 0.5 * phase_diff_shift, phase_diff_num);
x_pdf = zeros(phase_diff_num, 1);

for tr_id = 1:num_target_trajectories
    for p_id = 1 : size(xs, 1) - 1
		x = phase_diff(p_id, tr_id);
		x_id = floor((x - phase_diff_min) / (phase_diff_max - phase_diff_min + eps) * phase_diff_num) + 1;
		x_pdf(x_id) = x_pdf(x_id) + 1;
	end
end

x_pdf = x_pdf / (num_obs_periods * num_target_trajectories * phase_diff_shift);
norm = sum(x_pdf) * phase_diff_shift;
norm_diff = 1.0 - norm

fig = figure;
hLine = plot(phase_diff_bins, x_pdf);
set(gca, 'FontSize', 30);
xlabel('phase diff', 'Interpreter', 'latex');
xlim([phase_diff_bins(1) phase_diff_bins(end)])
set(gca, 'FontSize', 30);
ylabel('$PDF$', 'Interpreter', 'latex');

savefig(sprintf('%s/phase_diff_pdf_check_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/phase_diff_pdf_check_%s.pdf', home_figures_path, suffix_save));

dlmwrite(sprintf('%s/mean_trajectories_%s.txt', home_figures_path, suffix_save), mean_evo(:, 1:num_target_trajectories), 'delimiter', '\t')
dlmwrite(sprintf('%s/energy_trajectories_%s.txt', home_figures_path, suffix_save), energy_evo(:, 1:num_target_trajectories), 'delimiter', '\t')
dlmwrite(sprintf('%s/phase_diff_trajectories_%s.txt', home_figures_path, suffix_save), phase_diff, 'delimiter', '\t')


