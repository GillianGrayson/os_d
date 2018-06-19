clear all;

sys_id = 0;
task_id = 5;
prop_id = 0;

seed = 1;
mns = 1000000;
N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;
prm_E = 1.0;
prm_U = 0.15;
prm_J = 1.0;
start_type = 0;
start_state = 0;

T = 2 * pi / drv_freq;

cd_dump_deep = 1;
cd_num_sub_steps = 100;

size_sys = N + 1;

tr_id = 0;

is_mean = 1;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    seed, ...
    mns, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase, ...
    prm_E, ...
    prm_U, ...
    prm_J, ...
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
hLine = imagesc(dump_periods, states, adr);
title_str = sprintf('config(%d %d %d) tr(%d) N(%d) drv(%d %0.2f) prm(%0.2f %0.2f %0.2f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    tr_id, ...
    N, ... 
    drv_type, ...
    drv_ampl, ...
    prm_E, ...
    prm_U, ...
    prm_J);
title(title_str, 'FontSize', 33, 'interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([dump_periods(1) dump_periods(end)])
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$\rho_{n,n}$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
% set(gca, 'XAxisLocation', 'top');
% set(gca,'xticklabel',[]);
% set(gca, 'Position', [0.15 0.50 0.70 0.40])
hold all;

if(is_mean == 1)
    fn = sprintf('%s/mean_evo_%s.txt', data_path, suffix);
    mean_evo_data = importdata(fn);
    mean_evo = mean_evo_data(:, tr_id + 1);
    hLine = plot(dump_periods, mean_evo / size_sys);
end

propertyeditor(fig)


