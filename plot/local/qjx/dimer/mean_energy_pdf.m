clear all;

fig_path = 'E:/YandexDisk/Work/os_d/figures/dimer/qjx/sync_report';
data_path = '../../../../source/cpp/QJX/QJX';

sys_id = 0; 
task_id = 1;
prop_id = 1;
ss = 1;
mns = 1000000;
N = 100;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;
prm_E = 0.0;
prm_U = 0.15;
prm_J = 1.0;
start_type = 0;
start_state = 0;

n_tr = 100;

suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    ss, ...
    mns, ...
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
fn = sprintf('%s/energy_evo_%s.txt', data_path, suffix);
energy_evo = importdata(fn);

mean_all = [];
energy_all = [];
for tr_id = 1:n_tr
    mean_curr = mean_evo(:, tr_id);
    mean_all = vertcat(mean_all, mean_curr);
    energy_curr = energy_evo(:, tr_id);
    energy_all = vertcat(energy_all, energy_curr);
end

pdf2d.x_num_bins = 101;
pdf2d.y_num_bins = 101;
pdf2d.x_label = 'Mass center';
pdf2d.y_label = 'Energy';
pdf2d.x_bin_s = min(mean_all);
pdf2d.x_bin_f = max(mean_all);
pdf2d.y_bin_s = min(energy_all);
pdf2d.y_bin_f = max(energy_all);
pdf2d = oqs_pdf_2d_setup(pdf2d);
data2d = horzcat(mean_all, energy_all);
pdf2d = oqs_pdf_2d_update(pdf2d, data2d);
pdf2d = oqs_pdf_2d_release(pdf2d);
fig = oqs_pdf_2d_plot(pdf2d);
fn_fig = sprintf('%s/mean_energy_pdf_%s', fig_path, suffix);
oqs_save_fig(fig, fn_fig)