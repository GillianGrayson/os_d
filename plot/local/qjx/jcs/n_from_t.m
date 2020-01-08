clear all;

sys_id = 1;
task_id = 1;
prop_id = 0;
seed = 1;
mns = 1000000;

tr_id = 1;

N = 300;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
jcs_drv_part_1 = 0.5;
jcs_drv_part_2 = 0.5;
jcs_drv_ampl = 0.5;
jcs_prm_alpha = 5.0;
start_type = 0;
start_state = 0;

deep_dump = 0;
deep_num_steps = 50;


path = "../../../../source/cpp/QJX/QJX";

suffix = sprintf("setup(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)", ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    seed, ...
    mns, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    jcs_drv_part_1, ...
    jcs_drv_part_2, ...
    jcs_drv_ampl, ...
    jcs_prm_alpha, ...
    start_type, ...
    start_state);

fn = sprintf('%s/periods_%s.txt', path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

if(deep_dump == 1)
    num_periods = (num_dumps - 1) / (2 * deep_num_steps);
    dump_shift = 1 / (2 * deep_num_steps);
    dump_periods(1) = 0;
    for dump_id = 2:num_dumps
        dump_periods(dump_id) = dump_shift * (dump_id - 1);
    end
end

fn = sprintf('%s/mean_evo_%s.txt', path, suffix);
data = importdata(fn);

global_size = size(data, 1);

obs_evo = data(:, tr_id);

fig = figure;
hLine = plot(dump_periods, obs_evo, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$T/t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
hold all;


propertyeditor('on')