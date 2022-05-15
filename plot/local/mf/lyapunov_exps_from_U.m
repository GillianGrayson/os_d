clear all;

data_path = '../../../source/cpp/MF/MF';

task = 2;
U_start = 0.01;
U_shift = 0.01;
U_num = 50;
seed_start = 0;
seed_num = 1;
path = ''; 
mt = 1;
num_steps = 10000;
npt = 100;
np = 1000;
E = 0.0;
A = 3.4 ;
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
            seed-1);
        
        fn = sprintf('%s/exps_lpn_%s', data_path, fn_suffix);
        data = importdata(fn);
        
        for lpn_id = 1:num_lpns
            lpn_exps(U_id, lpn_id) = lpn_exps(U_id, lpn_id) + data(lpn_id) / seed_num;
        end
        
    end
    
end


fig = figure;

propertyeditor(fig);

for lpn_id = 1:num_lpns
    hLine = plot(Us, lpn_exps(:, lpn_id), 'LineWidth', 2);
    legend(hLine, sprintf('\lambda #%d', lpn_id));
    hold all;
end

set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
