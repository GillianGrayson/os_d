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

prm_E = 0.0;
prm_J = 1.0;

seed = 1;

U_start = 0.01;
U_shift = 0.01;
U_num = 100;

sys_size = N + 1;

Us = linspace(U_start, U_start + U_shift * (U_num - 1), U_num);
Ns = linspace(0, 1, sys_size);

rhonn_all = zeros(U_num, sys_size);

for U_id = 1:U_num
    
    U = Us(U_id)
    
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
    
    rhonn = diag(rho_mtx);
    norm_check = 1 - sum(rhonn)
    
    rhonn_all(U_id, :) = rhonn / max(rhonn);
end

fig = figure;
hLine = imagesc(Us, Ns, rhonn_all');
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n/N$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
h.Label.Interpreter = 'latex';
title(h, '\rho_{n,n}', 'FontSize', 33);
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

savefig(sprintf('%s/rhonn_from_U_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/rhonn_from_U_%s.pdf', home_figures_path, suffix_save));


