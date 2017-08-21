clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = '/qj_results/delta_0.1000';

data_path = sprintf('%s%s', data_path, prefix);

N = 500;

U_begin = 0.01;
U_step = 0.01;
U_num = 75;

num_seeds = 100;

Us = zeros(U_num, 1);
states = linspace(1, N, N) / N;

rho = zeros(U_num, N);

for U_id = 1:U_num
	
	U_id = U_id
	Us(U_id) = U_begin + U_step * (U_id - 1);
   
	curr_rho_avg = zeros(N,1);
   
	for seed = 1:num_seeds 
		path_to_folder = sprintf('%s/N_%d/U_%0.4f/rnd_%d', ...
			data_path, ...
			N, ...
			Us(U_id), ...
			seed-1);
       
      
		path = sprintf('%s/avg_diag_rho_trajectory_0.txt', path_to_folder);
		data = importdata(path);
       
		curr_rho_avg = curr_rho_avg + data;
	end
   
	curr_rho_avg = curr_rho_avg / num_seeds;
   
	norm = sum(curr_rho_avg)
   
	rho(U_id, :) = curr_rho_avg;
    
end

rho = rho';

fig = figure;
%hLine = contourf(Us, states, log10(rho + eps), 100, 'LineStyle', 'none');
hLine = imagesc(Us, states, rho + eps);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '<\rho_{n,n}>');
set(gca,'YDir','normal');

fn_suffix = sprintf('N(%d)', ...
        N);
		
savefig(sprintf('%s/bifurcation_%s.fig', home_figures_path, fn_suffix));