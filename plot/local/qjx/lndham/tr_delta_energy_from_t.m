clear all;

base_color = [0 0 1];
jump_color = [0 1 0];
renorm_color = [1 0 0];
linewidth = 2;

T = 1;

sys_id = 4;
task_id = 8;
prop_id = 0;

N = 10;
seed = 1;
alpha = 0.5;

lpn_type = -1;
delta_lpn_s = 1e-3;
delta_lpn_f_h = 1e-3;
delta_lpn_f_l = 1e-3;

deep_dump = 1;
deep_num_steps = 200;

num_target_traj = 10;
tr_id = 2;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(%d_%d_%d)_rnd(1_1000000)_N(%d)_rnd(%d)_alpha(%0.4f)_T(%0.4f)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    N, ...
    seed, ...
    alpha, ...
    T, ...
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

fn = sprintf('%s/lambda_evo_%s.txt', data_path, suffix);
lambda_evo_data = importdata(fn);

lambda_evo = lambda_evo_data(:, num_target_traj + tr_id);

spec_evo_base = spec_evo_data(:, 2 * (tr_id - 1) + 1);
spec_evo_var = spec_evo_data(:, 2 * (tr_id - 1 + num_target_traj) + 1);
spec_evo_diff = abs(spec_evo_base - spec_evo_var);

time_diff = dump_periods(2) - dump_periods(1);

fn = sprintf('%s/jump_times_%d_%s.txt', data_path, (tr_id - 1), suffix);
jump_data = importdata(fn);

interruptions = vertcat(jump_data, x_ticks);

fig = figure;
yyaxis left
propertyeditor(fig);
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$\Delta(t)$', 'Interpreter', 'latex');
set(gca, 'YScale', 'log')
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
hold all;
yyaxis right;
set(gca, 'FontSize', 30);
ylabel('$L(t) = \frac{1}{t} \sum \ln{\frac{\Delta_f}{\Delta_s}}$', 'Interpreter', 'latex');
h = plot(dump_periods, lambda_evo);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h = plot([dump_periods(1), dump_periods(end)], [0, 0], 'Color', 'k', 'LineStyle', '--');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


lambda_avg = mean(lambda_data(num_target_traj+1:end));
title(sprintf('$\\langle \\lambda \\rangle_{%d \\mathrm{traj}} = %0.2e \\quad \\lambda=%0.2e$', num_target_traj, lambda_avg, lambda_data(num_target_traj + tr_id)), 'Interpreter', 'latex')

    

renorm_diffs_sum_log = sum(log(spec_evo_diff(200:200:end) / delta_lpn_s))
