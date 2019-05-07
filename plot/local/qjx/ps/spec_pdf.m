clear all;

sys_id = 2;
task_id = 1;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 10;

diss_type = 1;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
ps_drv_part_1 = 0.98; 
ps_drv_part_2 = 1.00; 
ps_drv_ampl = 3.2;
ps_prm_alpha = 5;
ps_prm_d = 0.2500;
ps_prm_g = 0.2500;
ps_diss_w = 0.05;

start_type = 0;
start_state = 0;

num_bins = 201;

path = "../../../../source/cpp/QJX/QJX";

suffix = sprintf("rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
    seed, ...
    mns, ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    diss_type, ...
    ps_diss_w, ...
    ps_drv_part_1, ...
    ps_drv_part_2, ...
    ps_drv_ampl, ...
    ps_prm_alpha, ...
    ps_prm_d, ...
    ps_prm_g, ...
    start_type, ...
    start_state);

fn = sprintf('%s/spec_evo_%s.txt', path, suffix);
data = importdata(fn);

all_re = data(:, 1);
all_im = data(:, 2);
for tr_id = 1:num_trajectories-1
    all_re = vertcat(all_re, data(:, 1 + 2*tr_id));
    all_im = vertcat(all_im, data(:, 2 + 2*tr_id));
end
global_size = size(all_re, 1);

min_re = min(all_re);
max_re = max(all_re);
min_re = -3;
max_re = 3;
shift_re = (max_re - min_re) / num_bins;

min_im = min(all_im);
max_im = max(all_im);
min_im = -3;
max_im = 3;
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


