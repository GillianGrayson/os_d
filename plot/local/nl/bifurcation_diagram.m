clear all;

data_path = '../../../source/cpp/NL/NL';

N = 500;

U_start = 0.00;
U_shift = 0.005;
U_num = 150;

num_periods = 20000;

num_seeds = 10;

states = linspace(1, N, N);

Us = zeros(U_num, 1);
BD = zeros(U_num , N);

for U_id = 1 : U_num
    
    U = U_start + U_shift * (U_id - 1)
    Us(U_id) = U;
    
    for seed = 1:num_seeds
        
        fn_suffix = sprintf('U(%0.4f)_seed(%d).txt', ...
            U, ...
            seed-1);
        
        fn = sprintf('%s/data_%s', data_path, fn_suffix);
        data = importdata(fn);
        
        theta = data(:,1);
        phi = data(:,2);
        
        coordinate = N/2*(cos(theta)+1);
        
        for per_id = 1 : num_periods
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