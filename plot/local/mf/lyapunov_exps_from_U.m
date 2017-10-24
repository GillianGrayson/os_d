clear all;

data_path = '../../../source/cpp/MF/MF';

task = 0;
U_start = 0.03;
U_shift = 0.03;
U_num = 100;
seed_start = 0;
seed_num = 1;
path = ''; 
mt = 0;
num_steps = 1000;
npt = 2000;
np = 2000;
E = 1.0;
A = 1.5 ;
omega = 1.0;
phase = 0.0;
gamma = 0.1;
J = 1.0;

num_lpns = 2;

lpn_exps = zeros(U_num, num_lpns);

Us = zeros(U_num, 1);

for U_id = 1 : U_num
    
    U = U_start + U_shift * (U_id - 1)
    Us(U_id) = U;
    
    for seed = 1:seed_num
        
        fn_suffix = sprintf('mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed(%d).txt', ...
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
            lpn_exps(U_id, lpn_id) = lpn_exps(U_id, lpn_id) + data(end, lpn_id) / seed_num;
        end
        
    end
    
end


fig = figure;

propertyeditor(fig);

for lpn_id = 1:num_lpns
    hLine = plot(Us, lpn_exps(:, lpn_id), 'LineWidth', 2);
    legend(hLine, sprintf('\lambda #%d', lpn_id))
    hold all;
end

set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
