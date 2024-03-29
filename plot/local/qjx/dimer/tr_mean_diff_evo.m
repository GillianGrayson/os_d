clear all;

sys_id = 0;
task_id = 8;
prop_id = 0;

lpn = 1e-4;

tr_id = 5;
num_trajectories = 5;
num_plot_traj = 5;

seed = 1;
mns = 1000000;

N = 200;

diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;

drv_type = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;

prm_E = 1.0;
prm_U = 0.5;
prm_J = 1.0;

start_type = 0;
start_state = 0;

T = 2 * pi / drv_freq;

cd_dump_deep = 1;
cd_num_sub_steps = 100;

size_sys = N + 1;



data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(-1_%0.4f_%0.4f_%0.4f)', ...
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
    start_state, ...
    log10(lpn), ...
    log10(lpn), ...
    log10(lpn));

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

states = linspace(0, size_sys - 1, size_sys - 1) / (size_sys - 1);

fn = sprintf('%s/mean_evo_%s.txt', data_path, suffix);
mean_evo_data = importdata(fn);

fn = sprintf('%s/mean_lpn_evo_%s.txt', data_path, suffix);
mean_lpn_evo_data = importdata(fn);

fig = figure;

mean_evo_base = mean_evo_data(:, tr_id);
mean_evo_var = mean_evo_data(:, tr_id + num_trajectories);

mean_evo_diff = abs(mean_evo_base - mean_evo_var) / size_sys;

mean_evo_diff_periods = zeros(num_periods, 1);
p_id = 1;
for i = 1 : 2 * cd_num_sub_steps : size(mean_evo_diff, 1)
    mean_evo_diff_periods(p_id) = mean_evo_diff(i);
    p_id = p_id + 1;
end

pks = findpeaks(mean_evo_diff);

hLine = plot(dump_periods, mean_evo_diff);
legend(hLine, sprintf('tr=%d', tr_id))
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$|\Delta|$', 'Interpreter', 'latex');
hold all;
h = plot([dump_periods(1) dump_periods(end)],  [lpn lpn]);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(gca, 'YScale', 'log')
ylim([1e-12, 10]);


propertyeditor(fig)

