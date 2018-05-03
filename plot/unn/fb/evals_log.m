clear all;

drt = 0;
N = 501;
E = -1;
J = -1;

U = 0.5;

g = 0.1;
A = -1.5;
omega = 1;
seed = 1;

np = 1000;

start_id = 1;
end_id = N;

data_path = '../../../data/cluster/unn';

bin_begin = 1e-10;
num_decades = 10;
bin_end = bin_begin * 10.^num_decades;
num_bin_per_decade = 5;

num_bin = num_bin_per_decade * num_decades;

bin_borders = zeros(num_bin + 1, 1);
bin_centers = zeros(num_bin, 1);
bin_pdf = zeros(num_bin, 1);

for bin_id = 1 : num_bin + 1
    bin_borders(bin_id) = bin_begin * 10.^((bin_id - 1) / num_bin_per_decade);
    if (bin_id <= num_bin)
        bin_centers(bin_id) = bin_begin * 10.^((bin_id - 1 + 0.5) / num_bin_per_decade);
    end
end

bin_diff = diff(bin_borders);

non_inc_count = 0;

local_path = sprintf('np_%d/drt_%d/N_%d/E0_%0.4f/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_%d', ...
    np, ...
    drt, ...
    N-1, ...
    E, ...
    J, ...
    U, ...
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

curr_evals = curr_evals(start_id : end_id);

for eval_id = 1 : size(curr_evals, 1)
    
    curr_eval = curr_evals(eval_id);
    
    if ((curr_eval >= bin_begin) && (curr_eval <= bin_end))
        bin_id = floor((log10(curr_eval) - log10(bin_begin)) * num_bin / (log10(bin_end) - log10(bin_begin) + eps)) + 1;
        bin_pdf(bin_id) = bin_pdf(bin_id) + 1;
    else
        non_inc_count = non_inc_count + 1;
    end
    
end

non_inc_count = non_inc_count
norm = sum(bin_pdf)

for bin_id = 1 : num_bin
    bin_pdf(bin_id) = bin_pdf(bin_id) / (norm * bin_diff(bin_id));
end

norm_check = 0;
for bin_id = 1 : num_bin
    norm_check = norm_check + bin_pdf(bin_id) * bin_diff(bin_id);
end
norm_check = norm_check

fig = figure;
propertyeditor(fig);

hLine = plot(bin_centers, bin_pdf, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$\lambda$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$P(\lambda)$', 'Interpreter', 'latex');
set(gca,'XScale','log');
set(gca,'YScale','log');
propertyeditor('on');

