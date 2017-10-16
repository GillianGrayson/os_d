clear all;

data_path = '../../../source/cpp/MF/MF';

task = 0;
U_start = 0.03;
U_shift = 0.03;
U_num = 100;
seed_start = 0;
seed_num = 10;
path = ''; 
mt = 0;
num_steps = 1000;
npt = 2000;
np = 5000;
E = 1.0;
A = 1.5 ;
omega = 1.0;
phase = 0.0;
gamma = 0.1;
J = 1.0;

N = 1000;

states = linspace(1, N, N);

Us = zeros(U_num, 1);
BD = zeros(U_num , N);

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
        
        fn = sprintf('%s/nu_%s', data_path, fn_suffix);
        data = importdata(fn);
        
        nu = data(2:end,1);
        
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
    
end

BD = BD';

fig = figure;
hLine = imagesc(Us, states, BD);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
%title(h, '$PDF$', 'Interpreter', 'latex');
set(gca,'YDir','normal');