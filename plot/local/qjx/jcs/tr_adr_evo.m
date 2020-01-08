clear all;

tr_id = 0;

seed = 1;
num_seeds = 1000000;

N = 300;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;

T = 2.0;

drv_T_1 = 1.00 * T;
drv_T_2 = 1.00 * T;
drv_A = 3.0;

prm_alpha = 5.0;

start_type = 0;
start_state = 0;

T = drv_T_1 + drv_T_2;


cd_dump_deep = 1;
cd_num_sub_steps = 50;

is_mean = 1;


data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(1_4_0)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)', ...
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
    dump_periods = dump_periods / T;
end

adr = zeros(N, num_dumps);

fn = sprintf('%s/adr_%d_%s.txt', data_path, tr_id, suffix);
adr_data = importdata(fn);

for dump_id = 1:num_dumps
    for i = 1:N
        adr(i, dump_id) = adr_data((dump_id-1)*N + i);
    end
    
    curr_data = adr(:, dump_id);
    curr_max = max(curr_data);
    
    adr(:, dump_id) = curr_data / curr_max;
end

states = linspace(1, N, N) / N;

fig = figure;
hLine = imagesc(dump_periods, states, adr);
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$n/N$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$|\psi_n|^2$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
hold all;

if(is_mean == 1)
    fn = sprintf('%s/mean_evo_%s.txt', data_path, suffix);
    mean_evo_data = importdata(fn);
    mean_evo = mean_evo_data(:, tr_id + 1);
    hLine = plot(dump_periods, mean_evo / N);
end

propertyeditor(fig)


