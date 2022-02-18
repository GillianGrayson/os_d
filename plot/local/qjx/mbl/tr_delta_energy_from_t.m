clear all;

base_color = [0 0 1];
jump_color = [0 1 0];
renorm_color = [1 0 0];
linewidth = 2;

T = 10;

sys_id = 3;
task_id = 8;
prop_id = 0;

is_jump = 1;

seed = 1;
mns = 1000000;

Nc = 8;

W_seed = 10;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

prm_W = 20.0;
prm_U = 1.0;
prm_J = 1.0;

start_type = 0;
start_state = 0;

lpn_type = -1;
delta_lpn_s = 1e-6;
delta_lpn_f_h = 1e-6;
delta_lpn_f_l = 1e-6;

deep_dump = 1;
deep_num_steps = 200;

num_target_traj = 5;
tr_id = 1;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    seed, ...
    mns, ...
    Nc, ...
    W_seed, ...
    W_mns, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    prm_W, ...
    prm_U, ...
    prm_J, ...
    start_type, ...
    start_state, ...
    lpn_type, ...
    log10(delta_lpn_s), ...
    log10(delta_lpn_f_h), ...
    log10(delta_lpn_f_l));

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

num_periods = (size(dump_periods, 1) - 1) / deep_num_steps;
x_ticks = T * linspace(1, num_periods, num_periods)';
x_labels = {};
for i = 1 : num_periods
    x_labels{i} = sprintf('%d\\tau', i);
end

fn = sprintf('%s/spec_lpn_evo_%s.txt', data_path, suffix);
spec_evo_data = importdata(fn);

fn = sprintf('%s/lambda_%s.txt', data_path, suffix);
lambda_data = importdata(fn);


spec_evo_base = spec_evo_data(:, 2 * (tr_id - 1) + 1);
spec_evo_var = spec_evo_data(:, 2 * (tr_id - 1 + num_target_traj) + 1);
spec_evo_diff = abs(spec_evo_base - spec_evo_var);

time_diff = dump_periods(2) - dump_periods(1);

fn = sprintf('%s/jump_times_%d_%s.txt', data_path, (tr_id - 1), suffix);
jump_data = importdata(fn);

interruptions = vertcat(jump_data, x_ticks);

fig = figure;
propertyeditor(fig);
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$\Delta(t)$', 'Interpreter', 'latex');
set(gca, 'YScale', 'log')
%ylim([1e-12, 0.01]);
xtickangle(0)
xticks(x_ticks);
xticklabels(x_labels);
hold all;
curr_jump_id = 1;
prev_start = 1;
for p_id = 1:num_periods
    
    while(curr_jump_id <= size(jump_data, 1) && jump_data(curr_jump_id) < dump_periods(p_id * deep_num_steps))
        jump_time = jump_data(curr_jump_id);
        jump_time_id = floor(jump_time/time_diff);
        if(jump_time_id + 1 ~= (p_id-1) * deep_num_steps)
            plot_d = spec_evo_diff(prev_start:jump_time_id + 1);
            plot_t = dump_periods(prev_start:jump_time_id + 1);
            h = plot(plot_t, plot_d, 'LineWidth', linewidth, 'Color', base_color);
            flipud(h)
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            h = plot([jump_time jump_time], [spec_evo_diff(jump_time_id + 1) spec_evo_diff(jump_time_id + 2)], 'Color', jump_color, 'LineWidth', linewidth, 'LineStyle', ':');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
            ololo = 1;
        end

        prev_start = jump_time_id + 2;
        curr_jump_id = curr_jump_id + 1;
    end
    
    plot_d = spec_evo_diff(prev_start:p_id * deep_num_steps);
    plot_t = dump_periods(prev_start:p_id * deep_num_steps);
    h = plot(plot_t, plot_d, 'LineWidth', linewidth, 'Color', base_color);
    flipud(h)
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    if (p_id < num_periods)
        h = plot([dump_periods(p_id * deep_num_steps) + 0.5 * time_diff, dump_periods(p_id * deep_num_steps) + 0.5 * time_diff], ...
            [spec_evo_diff(p_id * deep_num_steps), spec_evo_diff(p_id * deep_num_steps + 1)], ...
            'LineWidth', linewidth, 'Color', renorm_color, 'LineStyle', ':');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        prev_start = p_id * deep_num_steps + 1;
    end  
end

h = plot([dump_periods(1), dump_periods(end)], [delta_lpn_s, delta_lpn_s], 'Color', 'k', 'LineStyle', '--');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
box on;

    

