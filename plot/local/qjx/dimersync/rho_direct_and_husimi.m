clear all;

fig_path = 'E:/YandexDisk/Work/os_d/figures/dimer/qjx/sync_report';
data_path = '../../../../source/cpp/QJX/QJX';

sys_id = 6; 
task_id = 1;
prop_id = 1;
ss = 1;
mns = 1000000;
N = 100;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
prm_E = 0.0;
prm_U = 0.125;
prm_J = 1.0;
start_type = 0;
start_state = 0;

drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;

drv_ampl_2 = 0.3;
drv_freq_2 = 0.6180339887498948;
drv_phase_2 = 0.0;

main_period = 1;

num_trajectories = 100;
num_periods = 10;
period_plot = num_periods;

sys_size = N + 1;

rho_total = zeros(sys_size);
count_rho = 0;

phi_size = 40;
phi_begin = 0;
phi_end = 2*pi;
phis = linspace(phi_begin, phi_end, phi_size)';

nu_size = 40;
nu_begin = 0;
nu_end = pi;
nus = linspace(nu_begin, nu_end, nu_size)';

for tr_id = 1:num_trajectories
    
    suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv1(%0.4f_%0.4f_%0.4f)_drv2(%0.4f_%0.4f_%0.4f)_T(%d)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        ss, ...
        mns, ...
        N, ...
        diss_type, ...
        diss_gamma, ...
        diss_phase, ...
        drv_ampl, ...
        drv_freq, ...
        drv_phase, ...
        drv_ampl_2, ...
        drv_freq_2, ...
        drv_phase_2, ...
        main_period, ...
        prm_E, ...
        prm_U, ...
        prm_J, ...
        start_type, ...
        start_state);
    
    fn = sprintf('%s/phi_evo_%d_%s.txt', data_path, tr_id-1, suffix);
    phi_evo_raw = importdata(fn);
    phi_evo_raw = phi_evo_raw(1 : (num_periods + 1) * sys_size, :);
    phi_evo = complex(phi_evo_raw(:, 1), phi_evo_raw(:, 2));
    
    curr_phi = phi_evo(period_plot * sys_size + 1 : (period_plot + 1) * sys_size);
    curr_norm = sum(abs(curr_phi).^2);
    curr_rho = zeros(sys_size);
    for i = 1:sys_size
        for j = 1:sys_size
            curr_rho(i,j) = dot(curr_phi(i), curr_phi(j)) / curr_norm;
        end
    end
    norm_check = 1 - trace(curr_rho)
    count_rho = count_rho + 1
    
    rho_total = rho_total + curr_rho;
    
end

rho_total = rho_total / count_rho;
norm_check_result = 1 - trace(curr_rho)

fig = figure;
propertyeditor(fig);
states = linspace(1, sys_size, sys_size)';
hLine = imagesc(states, states, abs(rho_total));
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$m$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
h.Label.Interpreter = 'latex';
title(h, '$|\rho_{n,m}|$', 'FontSize', 33, 'Interpreter', 'latex');
set(gca,'YDir','normal');
fn_fig = sprintf('%s/rho_direct_%d_%s', fig_path, period_plot, suffix);
oqs_save_fig(fig, fn_fig)

tic
hus = husimi(nus, phis, rho_total);
toc 

fig = figure;
propertyeditor(fig);
hLine = imagesc(nus, phis, real(hus'));
set(gca, 'FontSize', 30);
xlabel('$\vartheta$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\varphi$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'H', 'FontSize', 33, 'Interpreter', 'latex');
set(gca,'YDir','normal');
hold all;
fn_fig = sprintf('%s/rho_husimi_%d_%s', fig_path, period_plot, suffix);
oqs_save_fig(fig, fn_fig)