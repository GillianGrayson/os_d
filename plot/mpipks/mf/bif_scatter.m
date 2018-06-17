clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'mf_results';

data_path = sprintf('%s/%s', data_path, prefix);

task = 2;
U_start = 0.05;
U_shift = 0.05;
U_num = 20;
seed_start = 0;
seed_num = 10;
path = ''; 
mt = 0;
num_steps = 10000;
npt = 2000;
np = 2000;
E = 0.0;
A = 0.0 ;
omega = 1.0;
phase = 0.0;
gamma = 0.1;
J = -1.0;

nppp = 100;

N = 100;

states = linspace(1, N, N);

Us = zeros(U_num, 1);
BD = zeros(U_num , N);

scatter_data_x = [];
scatter_data_y = [];

for U_id = 1 : U_num
    
    U = U_start + U_shift * (U_id - 1)
    Us(U_id) = U;
    
    for seed = 1:seed_num
        
        path_to_folder = sprintf('%s/task_%d/np_%d/nppp_%d/mt_%d/omega_%0.4f/phase_%0.4f/g_%0.4f/J_%0.4f/E_%0.4f/A_%0.4f/U_%0.4f/eps_0.0000000100/m_4/seed_%d', ...
			data_path, ...
			task, ...
			np, ...
			nppp, ...
            mt, ...
			omega, ...
			phase, ...
			gamma, ...
			J, ...
			E, ...
			A, ...
			Us(U_id), ...
			seed);
        
        fn_suffix = sprintf('t(%d)_mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed(%d).txt', ...
            task, ...
			mt, ...
            omega, ...
            phase, ...
            gamma, ...
            J, ...
            E, ...
            A, ...
            U, ...
            seed);
        
        fn = sprintf('%s/data_%s', path_to_folder, fn_suffix);
        data = importdata(fn);
        
        nu = data(2:end, 1);
        
        coordinate = N/2*(cos(nu)+1);
        
        for per_id = 1 : np
            tmp = coordinate(per_id) / N * (N-1);
            id = floor(tmp) + 1;
            BD(U_id, id) = BD(U_id, id) + 1;
        end
        
        coordinate = 0;
        data = 0;
        
    end
    
    BD(U_id, :) = BD(U_id, :) / max(BD(U_id, :));
    
    [pks,locs] = findpeaks(BD(U_id, :));
    
    for loc_id = 1:size(locs, 2)
        scatter_data_x = vertcat(scatter_data_x, U);
        scatter_data_y = vertcat(scatter_data_y, locs(loc_id));
    end
end

fig = figure;
hLine = scatter(scatter_data_x, scatter_data_y);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');

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

savefig(sprintf('%s/mf_bif_scatter_%s.fig', home_figures_path, fn_suffix));
