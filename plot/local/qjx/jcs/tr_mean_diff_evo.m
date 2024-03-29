clear all;

seed = 1;
num_seeds = 1000000;

N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;

T_part = 2;

drv_T_1 = 0.98 * T_part;
drv_T_2 = 1.00 * T_part;
drv_A = 2.0;

prm_alpha = 5.0;

start_type = 0;
start_state = 0;

T = drv_T_1 + drv_T_2;

cd_dump_deep = 1;
cd_num_sub_steps = 100;

num_tr = 7;

lim = 0.3;

is_mean = 1;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)', ...
    seed, ...
    num_seeds, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    drv_T_1, ...
    drv_T_2, ...
    drv_A, ...
    prm_alpha, ...
    start_type, ...
    start_state);

fn = sprintf('%s/periods_%s.txt', data_path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

if(cd_dump_deep == 1)
    num_periods = (num_dumps - 1) / (2 * cd_num_sub_steps);
    dump_shift_1 = drv_T_1 / cd_num_sub_steps;
    dump_shift_2 = drv_T_2 / cd_num_sub_steps;
    curr_id = 1;
    dump_periods(curr_id) = 0;
    for period_id = 1:num_periods
      
        for dump_id = 1:cd_num_sub_steps
            dump_periods(curr_id + 1) = dump_periods(curr_id) + dump_shift_1;
            curr_id = curr_id + 1;
        end
        
        for dump_id = 1:cd_num_sub_steps
            dump_periods(curr_id + 1) = dump_periods(curr_id) + dump_shift_2;
            curr_id = curr_id + 1;
        end

    end
end

dump_periods = dump_periods / T;

states = linspace(1, N, N) / N;

fn = sprintf('%s/mean_evo_%s.txt', data_path, suffix);
mean_evo_data = importdata(fn) / N;

mean_evo_base = mean_evo_data(:, 1);


fig = figure;
for tr_id = 2:num_tr
    mean_evo_var = mean_evo_data(:, tr_id);

    mean_evo_diff = abs(mean_evo_base - mean_evo_var);

    hLine = plot(dump_periods, mean_evo_diff);
    legend(hLine, sprintf('tr=%d', tr_id))
    set(gca, 'FontSize', 30);
    xlabel('$t/T$', 'Interpreter', 'latex');
    xlim([dump_periods(1) dump_periods(end)])
    set(gca, 'FontSize', 30);
    ylabel('$|\Delta|$', 'Interpreter', 'latex');
    hold all;
    
end

hLine = plot([min(dump_periods) max(dump_periods)], [lim lim]);
hold all;

propertyeditor(fig)


