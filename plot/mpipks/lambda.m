clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = '/qj_results/';

data_path = sprintf('%s%s', data_path, prefix);

task = 4;
delta = 0.01;
tt = 1000;
E = 0.0;
T = 2*pi;
A = 0.0;
N = 100;

U_begin = 0.01;
U_step = 0.01;
U_num = 100;

J = -1.0;
g = 0.1;

lambda_id = 1000;

num_seeds = 100;

Us = zeros(U_num, 1);

lambdas = zeros(U_num, 1);

for U_id = 1:U_num
	
	U_id = U_id

	Us(U_id) = U_begin + U_step * (U_id - 1);
   
	curr_lambda_avg = 0;
   
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
			Us(U_id), ...
			J, ...
			g, ...
			seed-1);
			
		path = sprintf('%s/lambda_evo_trajectory_1.txt', path_to_folder);
		data = importdata(path);
       
		curr_lambda = data(lambda_id);
		seed = seed;
	   
		curr_lambda_avg = curr_lambda_avg + data(lambda_id);
		
	end
   
	curr_lambda_avg = curr_lambda_avg / num_seeds;
   
	lambdas(U_id) = curr_lambda_avg;
    
end

fig = figure;
hLine = plot(Us, lambdas, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');

fn_suffix = sprintf('task(%d)_delta(%0.4f)_tt(%d)_E(%0.4f)_T(%0.4f)_A(%0.4f)_N(%d)_J(%0.4f)_g(%0.4f)', ...
		task, ...
		delta, ...
		tt, ...
		E, ...
		T, ...
		A, ...
		N, ...
		J, ...
		g);
		
savefig(sprintf('%s/lambda_%s.fig', home_figures_path, fn_suffix));