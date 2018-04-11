clear all;

sys_id = 1;
task_id = 4;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 1;

N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
jcs_drv_part_1 = 0.98;
jcs_drv_part_2 = 1.0;
jcs_drv_ampl = 1.44;
jcs_prm_alpha = 5.0;
start_type = 0;
start_state = 0;

deep_dump = 1;
deep_num_steps = 128;


path = "../../../../source/cpp/QJX/QJX";

suffix = sprintf("config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)", ...
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

fn = sprintf('%s/spec_evo_%s.txt', path, suffix);
data = importdata(fn);

global_size = size(data, 1);

spec_re = data(:, 1);
spec_im = data(:, 2);

fig = figure;
hLine = plot(dump_periods, spec_re, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$T/t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$Re(\theta)$', 'Interpreter', 'latex');
hold all;


propertyeditor('on')