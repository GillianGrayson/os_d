clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = '/qj_results/delta_0.1000';

data_path = sprintf('%s%s', data_path, prefix);

N = 500;

lambda_id = 1000;

U_begin = 0.01;
U_step = 0.01;
U_num = 75;

num_seeds = 100;

Us = zeros(U_num, 1);

lambdas = zeros(U_num, 1);

for U_id = 1:U_num
	
	U_id = U_id

	Us(U_id) = U_begin + U_step * (U_id - 1);
   
	curr_lambda_avg = 0;
   
	for seed = 1:num_seeds 
		path_to_folder = sprintf('%s/N_%d/U_%0.4f/rnd_%d', ...
			data_path, ...
			N, ...
			Us(U_id), ...
			seed - 1);
      
		path = sprintf('%s/lambda_evo_trajectory_1.txt', path_to_folder);
		data = importdata(path);
       
		curr_lambda = data(lambda_id)
		seed = seed
	   
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

fn_suffix = sprintf('N(%d)', ...
        N);
		
savefig(sprintf('%s/lambda_%s.fig', home_figures_path, fn_suffix));