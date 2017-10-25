clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'mf_results';

data_path = sprintf('%s/%s', data_path, prefix);

task = 3;
U_start = 0.03;
U_shift = 0.03;
U_num = 100;
seed_start = 0;
seed_num = 50;
path = ''; 
mt = 0;
num_steps = 10000;
npt = 2000;
np = 2000;
E = 1.0;
A = 1.5 ;
omega = 1.0;
phase = 0.0;
gamma = 0.1;
J = 1.0;

num_lpns = 3;

lpn_exps = zeros(U_num, num_lpns);

Us = zeros(U_num, 1);

for U_id = 1 : U_num
    
    U = U_start + U_shift * (U_id - 1)
    Us(U_id) = U;
    
    for seed = 1:seed_num
        
        path_to_folder = sprintf('%s/task_%d/mt_%d/omega_%0.4f/phase_%0.4f/g_%0.4f/J_%0.4f/E_%0.4f/A_%0.4f/U_%0.4f/seed_%d', ...
			data_path, ...
			task, ...
            mt, ...
			omega, ...
			phase, ...
			gamma, ...
			J, ...
			E, ...
			A, ...
			Us(U_id), ...
			seed);
        
        fn_suffix = sprintf('mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed(%d).txt', ...
            mt, ...
            omega, ...
            phase, ...
            gamma, ...
            J, ...
            E, ...
            A, ...
            U, ...
            seed);
        		
		fn = sprintf('%s/exps_lpn_%s', path_to_folder, fn_suffix);
        data = importdata(fn);
        
        for lpn_id = 1:num_lpns
            lpn_exps(U_id, lpn_id) = lpn_exps(U_id, lpn_id) + data(lpn_id) / seed_num;
        end
        
    end
    
end

fig = figure;

for lpn_id = 1:num_lpns
    hLine = plot(Us, lpn_exps(:, lpn_id), 'LineWidth', 2);
    legend(hLine, sprintf('%d', lpn_id));
    hold all;
end

set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');

set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');

fn_suffix = sprintf('mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed_num(%d)', ...
    mt, ...
    omega, ...
    phase, ...
    gamma, ...
    J, ...
    E, ...
    A, ...
    U, ...
    seed_num);

savefig(sprintf('%s/mf_lpn_exps_%s.fig', home_figures_path, fn_suffix));
