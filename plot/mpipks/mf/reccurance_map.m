clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'mf_results';

data_path = sprintf('%s/%s', data_path, prefix);

task = 3;
U = 1.56;
seed_start = 0;
seed_num = 50;
path = ''; 
mt = 0;
num_steps = 10000;
npt = 10000;
np = 10000;
E = 1.0;
A = 1.5 ;
omega = 1.0;
phase = 0.0;
gamma = 0.1;
J = 1.0;

N = 1000;

num_int = N;
x_n_begin = 0;
x_n_end = N;
x_n_shift = (x_n_end - x_n_begin) / num_int;
x_n_int = zeros(num_int, 1);
for int_id = 1:num_int
    x_n_int(int_id) = x_n_begin + int_id * x_n_shift - 0.5 * x_n_shift;
end
eps = 1.0e-6;

x_n_pdf = zeros(num_int, num_int);

num_hits = 0;

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
        U, ...
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
    
    fn = sprintf('%s/data_%s', path_to_folder, fn_suffix);
    data = importdata(fn);
    
    nu = data(2:end, 2);
    
    coordinate = N/2*(cos(nu)+1);
    
    for period_id = 1 : (np - 1)
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

num_possible_hits = seed_num * (np - 1);
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

savefig(sprintf('%s/mf_rec_map_%s.fig', home_figures_path, fn_suffix));
