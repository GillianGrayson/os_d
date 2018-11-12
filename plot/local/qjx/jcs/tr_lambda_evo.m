clear all;

tr_id = 2;

seed = 1;
num_seeds = 1000000;

N = 50;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;

T_part = 2;

drv_T_1 = 0.98 * T_part;
drv_T_2 = 1.00 * T_part;
drv_A = 4.0;

prm_alpha = 5.0;

start_type = 0;
start_state = 0;

T = drv_T_1 + drv_T_2;

cd_dump_deep = 0;
cd_num_sub_steps = 10;

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
    dump_shift = 1 / (2 * cd_num_sub_steps);
    dump_periods(1) = 0;
    for dump_id = 2:num_dumps
        dump_periods(dump_id) = dump_shift * (dump_id - 1);
    end
end


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


