clear all;


home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = '/check';

data_path = sprintf('%s%s', data_path, prefix);

N = 1000;
U = 0.7;
local_path = sprintf('%s/N_%d/U_%0.4f', ...
    data_path, ...
    N, ...
    U);

num_periods = 1000;

rnd_start = 0;
num_rnd = 100;

num_int = 200;
x_n_begin = 0;
x_n_end = 1000;
x_n_shift = (x_n_end - x_n_begin) / num_int;
x_n_int = zeros(num_int, 1);
for int_id = 1:num_int
    x_n_int(int_id) = x_n_begin + int_id * x_n_shift - 0.5 * x_n_shift;
end
eps = 1.0e-6;

x_n_pdf = zeros(num_int, num_int);

num_existed = 0;
num_hits = 0;
for rnd_cur = rnd_start : rnd_start + (num_rnd - 1)

	rnd_cur = rnd_cur

    path_to_folder = sprintf('%s/rnd_%d', local_path, rnd_cur);
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
num_lost_hits = num_possible_hits - num_hits

x_n_pdf = x_n_pdf / (num_hits * x_n_shift * x_n_shift);
x_n_pdf = x_n_pdf';

fig = figure;
hLine = imagesc(x_n_int, x_n_int, x_n_pdf);
set(gca, 'FontSize', 30);
xlabel('$x_n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$x_{n+1}$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF');
set(gca,'YDir','normal');
set(gca, 'Position', [0.10 0.15 0.75 0.75])

fn_suffix = sprintf('N(%d)_U(%0.4f).fig', ...
        N, ...
        U);
		
savefig(sprintf('%s/recurrence_map_%s.fig', home_figures_path, fn_suffix));