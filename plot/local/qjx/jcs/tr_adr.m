clear all;

tr_id = 6;

seed = 1;
num_seeds = 1000000;

N = 200;
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

adr = zeros(size_sys, num_dumps);

fn = sprintf('%s/adr_%d_%s.txt', data_path, tr_id, suffix);
adr_data = importdata(fn);

for dump_id = 1:num_dumps
    for i = 1:size_sys
        adr(i, dump_id) = adr_data((dump_id-1)*size_sys + i);
    end
end

states = linspace(1, size_sys, size_sys) / size_sys;


fig = figure;
hLine = plot(states, adr(:, end));
title_str = sprintf('$tr=%d$ $N=%d$ $E=%0.2f$ $U=%0.2f$ $J=%0.2f$ ', tr_id, N, prm_E, prm_U, prm_J);
title(title_str, 'FontSize', 33, 'interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
xlim([states(1) states(end)])
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
hold all;

propertyeditor(fig)


