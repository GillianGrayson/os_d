clear all;

seed = 1;
mns = 1000000;

tr_id = 1;

Nc = 8;

W_seed = 10;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

prm_W = 2;
prm_U = 1.0;
prm_J = 1.0;

start_type = 0;
start_state = 49;

num_bins = 1001;

path = "../../../../source/cpp/QJX/QJX";

suffix = sprintf('rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
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
    start_state);

fn = sprintf('%s/imbalance_evo_%s.txt', path, suffix);
data_x = importdata(fn);

fn = sprintf('%s/spec_evo_%s.txt', path, suffix);
data_y = importdata(fn);

all_re = data_x(:, tr_id);
all_im = data_y(:, 2 * (tr_id - 1) + 1);

global_size = size(all_re, 1);

min_re = min(all_re);
max_re = max(all_re);
shift_re = (max_re - min_re) / num_bins;

min_im = min(all_im);
max_im = max(all_im);
shift_im = (max_im - min_im) / num_bins;

int_re = zeros(num_bins, 1);
int_im = zeros(num_bins, 1);
for int_id = 1:num_bins
    int_re(int_id) = min_re + shift_re * int_id - 0.5 * shift_re;
    int_im(int_id) = min_im + shift_im * int_id - 0.5 * shift_im;
end

pdf = zeros(num_bins, num_bins);
for d_id = 1:global_size
    curr_re = all_re(d_id, 1);
    curr_im = all_im(d_id, 1);
    
    int_re_id = floor((curr_re - min_re) * num_bins / (max_re - min_re + 10e-8)) + 1;
    int_im_id = floor((curr_im - min_im) * num_bins / (max_im - min_im + 10e-8)) + 1;
    
    pdf(int_re_id, int_im_id) = pdf(int_re_id, int_im_id) + 1;
end

pdf = pdf / (global_size * shift_re * shift_im);

norm = sum(sum(pdf)) * shift_re * shift_im

fig = figure;
hLine = imagesc(int_re, int_im, log10(pdf' + 1e-8));
set(gca, 'FontSize', 30);
xlabel('$I$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\varepsilon$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$\log_{10}PDF$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
set(gca, 'Position', [0.125 0.15 0.75 0.75]);
set(h, 'Position', [0.885 0.15 0.02 0.75]);
hold all;

propertyeditor('on')



