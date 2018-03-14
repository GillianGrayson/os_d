clear all;

sys_id = 0;
task_id = 4;
prop_id = 0;

seed = 1;
mns = 1000000;
N = 100;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;
prm_E = 1.0;
prm_U = 0.01;
prm_J = 1.0;
start_type = 0;
start_state = 0;

T = 2 * pi / drv_freq;

cd_dump_deep = 1;
cd_num_sub_steps = 128;

size_sys = N + 1;

tr_id = 0;

is_mean = 0;

data_path = '../../../source/cpp/QJX/QJX';

suffix = sprintf('config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    seed, ...
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

fn = sprintf('%s/jump_times_%d_%s.txt', data_path, tr_id, suffix);
jump_times = importdata(fn);

fn = sprintf('%s/jump_norms_%d_%s.txt', data_path, tr_id, suffix);
jump_norms = importdata(fn);

fn = sprintf('%s/jump_etas_%d_%s.txt', data_path, tr_id, suffix);
jump_etas = importdata(fn);

fig = figure;
hLine = scatter(jump_times / T, jump_norms);
title_str = sprintf('config(%d %d %d) tr(%d) N(%d) drv(%d %0.2f) prm(%0.2f %0.2f %0.2f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    tr_id, ...
    N, ... 
    drv_type, ...
    drv_ampl, ...
    prm_E, ...
    prm_U, ...
    prm_J);
title(title_str, 'FontSize', 33, 'interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$||\psi||^2$', 'Interpreter', 'latex');
hold all;
hLine = scatter(jump_times / T, jump_etas);

propertyeditor(fig)


