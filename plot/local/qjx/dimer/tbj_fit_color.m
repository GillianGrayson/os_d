clear all;

data_path = '../../../../data/cluster/mpipks';

sys_id = 0;
task_id = 1;
prop_id = 0;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 1337;
num_obs_periods = 10000;

N = 100;

diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;

U_begin = 0.02;
U_step = 0.02;
U_num = 100;
Us = linspace(U_begin, U_num * U_step, U_num);

A_begin = 0.02;
A_step = 0.02;
A_num = 100;
As = linspace(A_begin, A_num * A_step, A_num);

alphas = zeros(U_num, A_num);
decades = zeros(U_num, A_num);

for A_id = 1:A_num
    
    A = A_begin + A_step * (A_id - 1);
    
    fn = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/N_%d/diss_%d_%0.4f_%0.4f/A_%0.4f.txt', ...
        data_path, ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        ex_deep, ...
        rk_ns, ...
        num_tp_periods, ...
        num_obs_periods, ...
        N, ...
        diss_type, ...
        diss_gamma, ...
        diss_phase, ...
        A);
    
    data = importdata(fn);
    
    alphas(:, A_id) = data(:, 1);
    decades(:, A_id) = data(:, 2);
end

fig = figure;
hLine = imagesc(Us, As, alphas');
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$A$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\alpha');
set(gca,'YDir','normal');

fig = figure;
hLine = imagesc(Us, As, decades');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$T$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'decades');
set(gca,'YDir','normal');
