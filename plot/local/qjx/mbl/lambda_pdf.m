clear all;

num_traj = 200;

lambda_num_bins = 400;
lambda_min = -0.25;
lambda_max = 0.25;

sys_id = 3;
task_id = 7;
prop_id = 0;

seed = 1;
mns = 1000000;

Nc = 8;

W_seed = 10;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

prm_W = 20.0;
prm_U = 1.0;
prm_J = 1.0;

start_type = 0;
start_state = 49;

delta_lpn_s = 1e-6;
delta_lpn_f_h = 1e-3;
delta_lpn_f_l = 1e-12;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    seed, ...
    mns, ...
    Nc, ...
    W_seed, ...
    W_mns, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    prm_W, ...
    prm_U, ...
    prm_J, ...
    start_type, ...
    start_state, ...
    0, ...
    log10(delta_lpn_s), ...
    log10(delta_lpn_f_h), ...
    log10(delta_lpn_f_l));

fn = sprintf('%s/lambda_%s.txt', data_path, suffix);
lambda_data = importdata(fn);


curr_lambdas = lambda_data(num_traj / 2  + 1 : end);



lambda_shift = (lambda_max - lambda_min) / lambda_num_bins;
lambda_bins = linspace(lambda_min + 0.5 * lambda_shift, lambda_max - 0.5 * lambda_shift, lambda_num_bins);
lambdas_pdf_deep = zeros(lambda_num_bins, 1);

non_inc_deep = 0;
for l_id = 1: size(curr_lambdas, 1)
    curr_l = curr_lambdas(l_id);
    if curr_l < lambda_max && curr_l > lambda_min
        curr_id = floor((curr_l - lambda_min) / (lambda_max - lambda_min) * lambda_num_bins) + 1;
        lambdas_pdf_deep(curr_id) = lambdas_pdf_deep(curr_id) + 1;
    else
        non_inc_deep = non_inc_deep + 1;
    end
end
non_inc_deep = non_inc_deep

sum_lambda_pdf_deep = sum(lambdas_pdf_deep);
lambdas_pdf_deep = lambdas_pdf_deep / (sum_lambda_pdf_deep * lambda_shift);
norm = sum(lambdas_pdf_deep) * lambda_shift;
norm_diff_deep = 1.0 - norm

fig = figure;
hLine = plot(lambda_bins, lambdas_pdf_deep);
set(gca, 'FontSize', 30);
xlabel('$\lambda$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('PDF', 'Interpreter', 'latex');
hold all;