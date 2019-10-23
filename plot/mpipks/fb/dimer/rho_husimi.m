clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/fb/dimer';

task = 0;

num_steps = 10000;
num_periods_trans = 100;
num_periods_obser = 0;

N = 100;

diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;

drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;

U = 0.1;
prm_E = 0.0;
prm_J = 1.0;

seed = 1;

sys_size = N + 1;

phi_size = 100;
nu_size = 100;

phi_begin = 0;
phi_end = 2*pi;

nu_begin = 0;
nu_end = pi;

phis = linspace(phi_begin, phi_end, phi_size)';
nus = linspace(nu_begin, nu_end, nu_size)';

path_to_folder = sprintf('%s/task_%d/int_%d_%d_%d/N_%d/diss_%d_%0.4f_%0.4f/drv_%d_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/seed_%d', ...
    data_path, ...
    task, ...
    num_steps, ...
    num_periods_trans, ...
    num_periods_obser, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase, ...
    prm_E, ...
    U, ...
    prm_J, ...
    seed);

fn = sprintf('%s/rho.txt', path_to_folder);
rho_data = importdata(fn);

rho_mtx = zeros(sys_size);
for s_id = 1 : size(rho_data, 1)
    curr_row = rho_data(s_id, 1);
    curr_col = rho_data(s_id, 2);
    rho_mtx(curr_row, curr_col) = rho_data(s_id, 3) + sqrt(-1) * rho_data(s_id, 4);
end

tic
hus = husimi(nus, phis, rho_mtx);
toc 

fig = figure;
hLine = imagesc(nus, phis, real(hus'));
set(gca, 'FontSize', 30);
xlabel('$\vartheta$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\varphi$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'H', 'FontSize', 33);
set(gca,'YDir','normal');
hold all;

suffix_save = sprintf('int(%d_%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_seed(%d)', ...
    num_steps, ...
    num_periods_trans, ...
    num_periods_obser, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase, ...
    prm_E, ...
    U, ...
    prm_J, ...
    seed);

savefig(sprintf('%s/rho_husimi_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rho_husimi_%s.pdf', home_figures_path, suffix_save));


