clear all;

sys_id = 1;
task_id = 4;
prop_id = 0;

seed = 1;
mns = 1000000;

N = 200;

diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;

jcs_drv_part_1 = 0.98;
jcs_drv_part_2 = 1.0;
jcs_drv_ampl = 3.2;

jcs_prm_alpha = 5.0;

start_type = 0;
start_state = 0;

cd_num_sub_steps = 100;
num_obs_periods = 100;

times = zeros(2 * cd_num_sub_steps * num_obs_periods + 1, 1);
for p_id = 1:num_obs_periods
    
    for ss_id = 1:cd_num_sub_steps
        index = (p_id - 1) *  2 * cd_num_sub_steps + ss_id + 1;
        times(index) = times(index - 1) + jcs_drv_part_1 / cd_num_sub_steps;
    end
    
    for ss_id = 1:cd_num_sub_steps
        index = (p_id - 1) *  2 * cd_num_sub_steps + cd_num_sub_steps + ss_id + 1;
        times(index) = times(index - 1) + jcs_drv_part_2 / cd_num_sub_steps;
    end 
end
times = times * jcs_prm_alpha;

tr_id = 0;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)', ...
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

fn = sprintf('%s/norm_evo_%s.txt', data_path, suffix);
norm_data = importdata(fn);
norm_evo = norm_data(:, tr_id + 1);

fn = sprintf('%s/jump_times_%d_%s.txt', data_path, tr_id, suffix);
jump_times = importdata(fn);
jump_times_after = jump_times + 10e-10;

fn = sprintf('%s/jump_norms_%d_%s.txt', data_path, tr_id, suffix);
jump_norms = importdata(fn);
jump_norms_after = ones(size(jump_norms, 1), 1);

fn = sprintf('%s/jump_etas_%d_%s.txt', data_path, tr_id, suffix);
jump_etas = importdata(fn);

all_times = vertcat(times, jump_times, jump_times_after);
all_norms = vertcat(norm_evo, jump_norms, jump_norms_after);

[all_times_sorted, order] = sort(all_times);
all_norms_sorted = zeros(size(all_times_sorted, 1), 1);

for t_id = 1:size(all_times_sorted, 1)
    all_norms_sorted(t_id) = all_norms(order(t_id));
end

fig = figure;
hLine = plot(all_times_sorted, all_norms_sorted);
title_str = sprintf('config(%d %d %d) tr(%d) N(%d) drv(%0.4f %0.4f %0.4f) prm(%0.2f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    tr_id, ...
    N, ... 
    jcs_drv_part_1, ...
    jcs_drv_part_2, ...
    jcs_drv_ampl, ...
    jcs_prm_alpha);
title(title_str, 'FontSize', 33, 'interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$||\psi||^2$', 'Interpreter', 'latex');
hold all;
hLine = scatter(jump_times, jump_norms, 'Marker', 'o', 'MarkerFaceColor', 'w');
set(get(get(hLine, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
hold all;
propertyeditor(fig)



