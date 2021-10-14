clear all;

base_color = [0 0 1];
jump_color = [0 1 0];
renorm_color = [1 0 0];
linewidth = 2;

T = 1;

sys_id = 5;
task_id = 8;
prop_id = 0;

N = 7;
seed = 1;
tau = 1;
k = -1;

lpn_type = -1;
delta_lpn = 1e-6;

deep_dump = 1;
deep_num_steps = 200;

num_target_traj = 50;
tr_id = 0;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(%d_%d_%d)_rnd(0_1000000)_N(%d)_seed(%d)_tau(%d)_k(%d)_T(%0.4f)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    N, ...
    seed, ...
    tau, ...
    k, ...
    T, ...
    lpn_type, ...
    log10(delta_lpn), ...
    log10(delta_lpn), ...
    log10(delta_lpn));

fn = sprintf('%s/periods_%s.txt', data_path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

if(deep_dump == 1)
    num_periods = (num_dumps - 1) / (1 * deep_num_steps);
    dump_shift = 1 / (1 * deep_num_steps);
    dump_periods(1) = 0;
    for dump_id = 2:num_dumps
        dump_periods(dump_id) = dump_shift * (dump_id - 1);
    end
end

dump_periods = dump_periods * T;

fn = sprintf('%s/norm_evo_%s.txt', data_path, suffix);
norm_evo_data = importdata(fn);
norm_evo = norm_evo_data(:, tr_id + 1);

fn = sprintf('%s/jump_times_%d_%s.txt', data_path, tr_id, suffix);
jump_times = importdata(fn);
jump_times = jump_times / T;
jump_times_after = jump_times + 10e-10;

fn = sprintf('%s/jump_norms_%d_%s.txt', data_path, tr_id, suffix);
jump_norms = importdata(fn);
jump_norms_after = ones(size(jump_norms, 1), 1);

fn = sprintf('%s/jump_etas_%d_%s.txt', data_path, tr_id, suffix);
jump_etas = importdata(fn);

all_times = vertcat(dump_periods, jump_times, jump_times_after);
all_norms = vertcat(norm_evo, jump_norms, jump_norms_after);

[all_times_sorted, order] = sort(all_times);
all_norms_sorted = zeros(size(all_times_sorted, 1), 1);

for t_id = 1:size(all_times_sorted, 1)
    all_norms_sorted(t_id) = all_norms(order(t_id));
end

fig = figure;
hLine = plot(all_times_sorted, all_norms_sorted);
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$||\psi||^2$', 'Interpreter', 'latex');
hold all;
hLine = scatter(jump_times, jump_norms);

propertyeditor(fig)


