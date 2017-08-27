clear all;

data_path = '../../../source/cpp/NL/NL';

N = 500;

U = 0.55;

num_periods = 2000;

seed_begin = 0;
seed_num = 100;

num_int = 500;
x_n_begin = 0;
x_n_end = 500;
x_n_shift = (x_n_end - x_n_begin) / num_int;
x_n_int = zeros(num_int, 1);
for int_id = 1:num_int
    x_n_int(int_id) = x_n_begin + int_id * x_n_shift - 0.5 * x_n_shift;
end
eps = 1.0e-6;

x_n_pdf = zeros(num_int, num_int);

num_hits = 0;
for seed = seed_begin : seed_begin + (seed_num - 1)
    
    fn_suffix = sprintf('U(%0.4f)_seed(%d).txt', ...
        U, ...
        seed);
    
    fn = sprintf('%s/data_%s', data_path, fn_suffix);
    data = importdata(fn);
    
    theta = data(:,1);
    phi = data(:,2);
    
    coordinate = N/2*(cos(theta)+1);
    
    for period_id = 1 : (num_periods - 1)
        x_n = coordinate(period_id);
        y_n = coordinate(period_id + 1);
        
        if x_n >= x_n_begin && x_n <= x_n_end && y_n >= x_n_begin && y_n <= x_n_end
            
            num_hits = num_hits + 1;
            
            x_id = floor((x_n - x_n_begin) * num_int / (x_n_end - x_n_begin + eps)) + 1;
            y_id = floor((y_n - x_n_begin) * num_int / (x_n_end - x_n_begin + eps)) + 1;
            
            x_n_pdf(x_id, y_id) = x_n_pdf(x_id, y_id) + 1;
        end
    end
end

num_possible_hits = seed_num * (num_periods - 1);
num_lost_hits = num_possible_hits - num_hits

x_n_pdf = x_n_pdf / (num_hits * x_n_shift * x_n_shift);
x_n_pdf = x_n_pdf';

fig = figure;
hLine = imagesc(x_n_int, x_n_int, log10(x_n_pdf + 1.0e-8));
set(gca, 'FontSize', 30);
xlabel('$x_n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$x_{n+1}$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$\log_{10}PDF$', 'Interpreter', 'latex');
set(gca,'YDir','normal');

propertyeditor('on')
