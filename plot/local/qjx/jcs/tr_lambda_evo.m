clear all;

tr_id = 11;

lpn_type = -1;
lpn_delta_s = log10(1.0e-3);
lpn_delta_f_high = log10(1.0e-1);
lpn_delta_f_low = log10(1.0e-5);

seed = 1;
num_seeds = 1000000;

N = 300;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;

T_part = 2.00;

drv_T_1 = 1.00 * T_part;
drv_T_2 = 1.00 * T_part;
drv_A = 4.0;

prm_alpha = 5.0;

start_type = 0;
start_state = 0;

T = drv_T_1 + drv_T_2;

cd_dump_deep = 1;
cd_num_sub_steps = 50;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(1_8_0)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
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
    start_state, ...
    lpn_type, ...
    lpn_delta_s, ...
    lpn_delta_f_high, ...
    lpn_delta_f_low);

fn = sprintf('%s/periods_%s.txt', data_path, suffix);
dump_periods = importdata(fn) - 2;
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


fn = sprintf('%s/lambda_evo_%s.txt', data_path, suffix);
evo_data = importdata(fn);


evo_curr = evo_data(:, tr_id + 1);

hLine = plot(dump_periods, evo_curr);
legend(hLine, sprintf('tr=%d', tr_id))
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');
hold all;

propertyeditor('on')


