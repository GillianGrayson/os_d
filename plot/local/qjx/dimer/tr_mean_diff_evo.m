clear all;

sys_id = 0;
task_id = 5;
prop_id = 0;

seed = 1;
mns = 1000000;
N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;
prm_E = 1.0;
prm_U = 0.05;
prm_J = 1.0;
start_type = 0;
start_state = 0;

T = 2 * pi / drv_freq;

cd_dump_deep = 1;
cd_num_sub_steps = 100;

size_sys = N + 1;

num_tr = 1;

is_mean = 1;

data_path = '../../../../source/cpp/QJX/QJX';
%data_path = 'C:/Users/user/Desktop/New folder/U(0.05)/T1_ds(1e-6)'

suffix = sprintf('rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
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

states = linspace(1, size_sys, size_sys) / size_sys;

fn = sprintf('%s/mean_evo_%s.txt', data_path, suffix);
mean_evo_data = importdata(fn);

mean_evo_base = mean_evo_data(:, 1);


fig = figure;
for tr_id = 1:num_tr
    mean_evo_var = mean_evo_data(:, tr_id + 1);

    mean_evo_diff = abs(mean_evo_base - mean_evo_var) / size_sys;
    
    dump_periods = dump_periods(1: 1 + cd_num_sub_steps * 200);
    mean_evo_diff = mean_evo_diff(1: 1 + cd_num_sub_steps * 200);
    
    full_periods_x = dump_periods(2:2 * cd_num_sub_steps:size(dump_periods, 1));
    full_periods_y = mean_evo_diff(2:2 * cd_num_sub_steps:size(dump_periods, 1));

    hLine = plot(dump_periods, mean_evo_diff);
    legend(hLine, sprintf('tr=%d', tr_id))
    set(gca, 'FontSize', 30);
    xlabel('$t/T$', 'Interpreter', 'latex');
    xlim([dump_periods(1) dump_periods(end)])
    set(gca, 'FontSize', 30);
    ylabel('$|\Delta|$', 'Interpreter', 'latex');
    hold all;
    hLine = plot(full_periods_x, full_periods_y, 'o');
    
end

propertyeditor(fig)


