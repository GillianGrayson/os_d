clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'mf_results';

data_path = sprintf('%s/%s', data_path, prefix);

task = 1;

U_start = 0.01;
U_shift = 0.01;
U_num = 200;

A_start = 0.01;
A_shift = 0.01;
A_num = 200;

ss_start = 0;
ss_num = 1;

seed_start = 0;
seed_num = 10;

path = '';
mt = 0;
num_steps = 10000;
npt = 2100;
np = 2100;
E = 1.0;
J = 1.0;
omega = 1.0;
phase = 0.0;
gamma = 0.1;

nppp = 100;

Us = linspace(U_start, U_shift * U_num, U_num);
As = linspace(A_start, A_shift * A_num, A_num);

lpns = zeros(U_num, A_num);

for U_id = 1 : U_num
    
    U = U_start + U_shift * (U_id - 1)
    
    for A_id = 1 : A_num
        
        A = A_start + A_shift * (A_id - 1);
        
        lpn_avg = 0;
        
        for ss = 1:ss_num
            
            for curr_seed = 1:seed_num
                
                seed = ss_start + (ss - 1) * seed_num + curr_seed;
                
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
                    ss);
                
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
                
                fn = sprintf('%s/exps_lpn_%s', path_to_folder, fn_suffix);
                data = importdata(fn);
                
                lpn_avg = lpn_avg + data(1);
                
            end
        end
        
        lpn_avg = lpn_avg / (ss_num * seed_num);
        
        lpns(U_id, A_id) = lpn_avg;
        
    end
end


fig = figure;
hLine = imagesc(Us, As, lpns');
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$A$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$\lambda$', 'Interpreter', 'latex');
set(gca,'YDir','normal');

fn_suffix = sprintf('mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_sm(%d)', ...
    mt, ...
    omega, ...
    phase, ...
    gamma, ...
    J, ...
    E, ...
    ss_num * seed_num);

savefig(sprintf('%s/lpns_space_%s.fig', home_figures_path, fn_suffix));
