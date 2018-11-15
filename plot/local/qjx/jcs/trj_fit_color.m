clear all;

data_path = '../../../../data/cluster/mpipks';

sys_id = 1;
task_id = 1;
prop_id = 0;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 1337;
num_obs_periods = 10000;

N = 200;

diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;

ampl_begin = 0.025;
ampl_step = 0.025;
ampl_num = 200;
ampls = linspace(ampl_begin, ampl_num * ampl_step, ampl_num);

T_begin = 0.025;
T_step = 0.025;
T_num = 200;
Ts = linspace(T_begin, T_num * T_step, T_num);

alphas = zeros(ampl_num, T_num);
decades = zeros(ampl_num, T_num);

for T_id = 1:T_num
    
    T = T_begin + T_step * (T_id - 1);
    
    fn = sprintf('%s/T_%0.4f.txt', ...
        data_path, ...
        T);
    
    data = importdata(fn);
    
    alphas(:, T_id) = data(:, 1);
    decades(:, T_id) = data(:, 2);
end

fig = figure;
hLine = imagesc(ampls, Ts, alphas');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$T$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\alpha');
set(gca,'YDir','normal');

fig = figure;
hLine = imagesc(ampls, Ts, decades');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$T$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'decades');
set(gca,'YDir','normal');
