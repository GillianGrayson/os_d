clear all;

sys_id = 0;
task_id = 5;
prop_id = 1;

seed = 0;
mns = 1000000;
N = 100;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 1;
drv_ampl = 3.4;
drv_freq = 1.0;
drv_phase = 0.0;
prm_E = 0.0;
prm_U = 0.5;
prm_J = 1.0;
start_type = 0;
start_state = 0;

T = 2 * pi / drv_freq;

cd_dump_deep = 1;
cd_num_sub_steps = 128;

size_sys = N + 1;

tr_id = 1;

is_mean = 1;

data_path = '../../../../source/cpp/QJX/QJX';

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

fn = sprintf('%s/periods_%s.txt', data_path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

if(cd_dump_deep == 1)
    num_periods = (num_dumps - 1) / (2 * cd_num_sub_steps);
    dump_shift = 1 / (2 * cd_num_sub_steps);
    dump_periods(1) = 0;
    for dump_id = 2:num_dumps
        dump_periods(dump_id) = dump_shift * (dump_id - 1);
    end
end

adr = zeros(size_sys, num_dumps);

states = linspace(1, size_sys, size_sys) / size_sys;

fn = sprintf('%s/mean_evo_%s.txt', data_path, suffix);
mean_evo_data = importdata(fn);
mean_evo = mean_evo_data(:, tr_id + 1);

fig = figure;
hLine = plot(dump_periods, mean_evo / size_sys);
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
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
hold all;

propertyeditor(fig)


