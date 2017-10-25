clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = '/qj_results/';

data_path = sprintf('%s%s', data_path, prefix);

task = 3;
delta = 0.1;
tt = 1000;
E = 0.0;
T = 2*pi;
A = 0.0;
N = 50;

U = 1.0;

J = -1.0;
g = 0.1;

num_periods = 1000;

num_seeds = 100;

num_int = 50;
x_n_begin = 0;
x_n_end = 50;
x_n_shift = (x_n_end - x_n_begin) / num_int;
x_n_int = zeros(num_int, 1);
for int_id = 1:num_int
    x_n_int(int_id) = x_n_begin + int_id * x_n_shift - 0.5 * x_n_shift;
end
eps = 1.0e-6;

x_n_pdf = zeros(num_int, num_int);

num_existed = 0;
num_hits = 0;

for seed = 1:num_seeds
    		
	path_to_folder = sprintf('%s/task_%d/delta_%0.4f/tt_%d/E_%0.4f/T_%0.4f/A_%0.4f/N_%d/U_%0.4f/J_%0.4f/g_%0.4f/rnd_%d', ...
			data_path, ...
			task, ...
			delta, ...
			tt, ...
			E, ...
			T, ...
			A, ...
			N, ...
			U, ...
			J, ...
			g, ...
			seed-1);
    
    
    path = sprintf('%s/mean_evo_trajectory_0.txt', path_to_folder);
    
    if (exist(path, 'file') == 2)
        
		num_existed = num_existed + 1;
        
        data = importdata(path);
        
        for period_id = 1 : (num_periods - 1)
            x_n = data(period_id);
            y_n = data(period_id + 1);
            
            if x_n >= x_n_begin && x_n <= x_n_end && y_n >= x_n_begin && y_n <= x_n_end
                
                num_hits = num_hits + 1;
                
                x_id = floor((x_n - x_n_begin) * num_int / (x_n_end - x_n_begin + eps)) + 1;
                y_id = floor((y_n - x_n_begin) * num_int / (x_n_end - x_n_begin + eps)) + 1;
            
                x_n_pdf(x_id, y_id) = x_n_pdf(x_id, y_id) + 1; 
            end
        end
        
    end
    
end

num_possible_hits = num_existed * (num_periods - 1);
num_hits = num_hits
num_lost_hits = num_possible_hits - num_hits

x_n_pdf = x_n_pdf / (num_hits * x_n_shift * x_n_shift);
x_n_pdf = x_n_pdf';

fig = figure;
hLine = imagesc(x_n_int, x_n_int, log10(x_n_pdf + 1.0e-8));
set(gca, 'FontSize', 30);
xlabel('$x_n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$x_{n+1}$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$\log_{10}PDF$', 'Interpreter', 'latex');
set(gca,'YDir','normal');

fn_suffix = sprintf('task(%d)_delta(%0.4f)_tt(%d)_E(%0.4f)_T(%0.4f)_A(%0.4f)_N(%d)_U(%0.4f)_J(%0.4f)_g(%0.4f)', ...
		task, ...
		delta, ...
		tt, ...
		E, ...
		T, ...
		A, ...
		N, ...
		U, ...
		J, ...
		g);
		
savefig(sprintf('%s/rm_%s.fig', home_figures_path, fn_suffix));
