clear all;

drt = 0;
N = 501;
E = -1;
J = -1;

U_start = 0.005;
U_shift = 0.005;
U_num = 150;

num_first_skips = floor(N/3);
num_last_skips = floor(N/3);

g = 0.1;
A = -1.5;
omega = 1;
seed = 1;

np = 1000;

data_path = '../../../data/cluster/unn';

Us = zeros(U_num, 1);
Ns = linspace(0, 1, N);

ratios = zeros(U_num, N);

for U_id = 1:U_num
    
    U_id = U_id
    
    Us(U_id) = U_start + (U_id-1) * U_shift;
    
   local_path = sprintf('np_%d/drt_%d/N_%d/E0_%0.4f/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_%d', ...
        np, ...
        drt, ...
        N-1, ...
        E, ...
        J, ...
        Us(U_id), ...
        g, ...
        A, ...
        omega, ...
        seed);
    
    fn = sprintf('%s/%s/rho.txt', data_path, local_path);
    mtx_data = importdata(fn);
    mtx = zeros(N);
    
    for s_id = 1 : size(mtx_data, 1)
        curr_row = mtx_data(s_id, 1);
        curr_col = mtx_data(s_id, 2);
        mtx(curr_row, curr_col) = mtx_data(s_id, 3) + sqrt(-1) * mtx_data(s_id, 4);
    end
    
    curr_evals = eig(mtx);
    
    curr_evals = abs(curr_evals);
    
    curr_evals = sort(curr_evals);
    
    cutted_evals = curr_evals(num_first_skips + 1:size(curr_evals, 1) - num_last_skips);
    size_cutted = size(cutted_evals, 1);
    
    s_n = zeros(size_cutted-1, 1);
    r_n = zeros(size_cutted-2, 1);
    
    for s_id = 1 : (size_cutted-1)
        s_n(s_id) = cutted_evals(s_id + 1) - cutted_evals(s_id);
    end
    
    for s_id = 1 : (size_cutted-2)
        r_n(s_id) = s_n(s_id + 1) / s_n(s_id);
    end
    
    for s_id = 1 : (size_cutted-2)
        r_n(s_id) = min(r_n(s_id), 1.0/r_n(s_id));
    end
    
    r_n_avg = 0.0;
    for s_id = 1 : (size_cutted-2)
        r_n_avg = r_n_avg + r_n(s_id);
    end
    
    r_n_avg = r_n_avg / (size_cutted-2);
    
    ratios(U_id) = r_n_avg;
    
end

fig = figure;
propertyeditor(fig);

hLine = plot(Us, ratios);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$ratio$', 'Interpreter', 'latex');
hold all;

