clear all;


home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = '/check';

data_path = sprintf('%s%s', data_path, prefix);

N = 200;

U_begin = 0.01;
U_step = 0.01;
U_num = 99;

Us = zeros(U_num, 1);
states = linspace(1, N, N) / N;

rho = zeros(U_num, N);

rnd_cur = 0;

for U_id = 1:U_num
   U(U_id) = U_begin + U_step * (U_id - 1);
   
   local_path = sprintf('%s/N_%d/U_%0.4f', ...
        data_path, ...
        N, ...
        U(U_id));
    
    path_to_folder = sprintf('%s/rnd_%d', local_path, rnd_cur);
    
    path = sprintf('%s/avg_diag_rho_trajectory_0.txt', path_to_folder);
    data = importdata(path);
    
    rho(U_id, :) = data;
    
end

rho = rho';

fig = figure;
hLine = imagesc(Us, states, rho);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '<\rho_{n,n}>');
set(gca,'YDir','normal');
set(gca, 'Position', [0.10 0.15 0.75 0.75])


		
savefig(sprintf('%s/check.fig', home_figures_path));