clear all;

sys_id = 1;
task_id = 1;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 1;

N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
jcs_drv_part_1 = 0.98;
jcs_drv_part_2 = 1.0;
jcs_drv_ampl = 3.2;
jcs_prm_alpha = 5.0;
start_type = 0;
start_state = 0;

num_bins = 201;

path = "../../../../source/cpp/QJX/QJX";

suffix = sprintf("config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)", ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    seed, ...
    mns, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    jcs_drv_part_1, ...
    jcs_drv_part_2, ...
    jcs_drv_ampl, ...
    jcs_prm_alpha, ...
    start_type, ...
    start_state);

fn = sprintf('%s/spec_evo_%s.txt', path, suffix);
data = importdata(fn);

global_size = size(data, 1);

min_re = min(data(:, 1));
max_re = max(data(:, 1));
shift_re = (max_re - min_re) / num_bins;

min_im = min(data(:, 2));
max_im = max(data(:, 2));
shift_im = (max_im - min_im) / num_bins;

int_re = zeros(num_bins, 1);
int_im = zeros(num_bins, 1);
for int_id = 1:num_bins
    int_re(int_id) = min_re + shift_re * int_id - 0.5 * shift_re;
    int_im(int_id) = min_im + shift_im * int_id - 0.5 * shift_im;
end

pdf = zeros(num_bins, num_bins);
for d_id = 1:global_size
    curr_re = data(d_id, 1);
    curr_im = data(d_id, 2);
    
    int_re_id = floor((curr_re - min_re) * num_bins / (max_re - min_re + 10e-8)) + 1;
    int_im_id = floor((curr_im - min_im) * num_bins / (max_im - min_im + 10e-8)) + 1;
    
    pdf(int_re_id, int_im_id) = pdf(int_re_id, int_im_id) + 1;
end

pdf = pdf / (global_size * shift_re * shift_im);

norm = sum(sum(pdf)) * shift_re * shift_im

fig = figure;
hLine = imagesc(int_re, int_im, pdf');
set(gca, 'FontSize', 30);
xlabel('$Re(\theta)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$Im(\theta)$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
set(gca, 'Position', [0.125 0.15 0.75 0.75]);
set(h, 'Position', [0.885 0.15 0.02 0.75]);
hold all;


propertyeditor('on')


