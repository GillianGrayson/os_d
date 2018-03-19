clear all;

N = 101;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/diag_rho_deep_evo.txt', data_path);
rho_data = importdata(fn);

times_num = 200;
num_T = 10;

rho_evo = zeros(N, times_num);

states = linspace(1, N, N);
times = zeros(times_num * num_T, 1);

for t_id = 1:times_num * num_T
    times(t_id) = t_id / times_num;
    start = 1 + (t_id - 1) * N;
    rho_evo(:, t_id) = rho_data(start : start+N-1);
end



fig = figure;
hLine = imagesc(times, states, rho_evo);
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
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

propertyeditor(fig)


