clear all;

data_path = '../../../../data/cluster/unn/qjx';

sys_id = 2;
task_id = 1;
prop_id = 0;
seed = 1;
mns = 1000000;
num_trajectories = 32;
num_tp_periods = 1000;
num_obs_periods = 1000;
ex_deep = 16;
rk_ns = 10000;

num_runs = 20;

diss_type = 1;
ps_diss_w = 0.05;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
ps_drv_part_1 = 0.98; 
ps_drv_part_2 = 1.00; 
ps_drv_ampl = 3.2;
ps_prm_alpha = 5;
ps_prm_d = 0.1;
ps_prm_g = 0.1;


start_type = 0;
start_state = 0;


all_re = zeros(num_runs * num_trajectories * (num_obs_periods + 1), 1);
all_im = zeros(num_runs * num_trajectories * (num_obs_periods + 1), 1);

for run_id = 1:num_runs
    
    ss = (run_id - 1) * num_trajectories;
    
    path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/N_%d_%d/diss_%d_%0.4f/drv_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
        data_path, ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        ex_deep, ...
        rk_ns, ...
        num_tp_periods, ...
        num_obs_periods, ...
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
        start_state, ...
        ss);
    
    suffix = sprintf("rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
        ss, ...
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
    
    path = sprintf('%s/spec_evo_%s.txt', path_to_folder, suffix);
    data = importdata(path);
    
    for tr_id = 1:num_trajectories
        start_run_id = (run_id-1) * num_trajectories * (num_obs_periods + 1) + 1;
        start_traj_id = start_run_id + (tr_id - 1) * (num_obs_periods + 1);
        end_traj_id = start_traj_id + (num_obs_periods + 1) - 1;
        all_re(start_traj_id:end_traj_id) = data(:, -1 + 2*tr_id);
        all_im(start_traj_id:end_traj_id) = data(:, 0 + 2*tr_id);
    end
end

num_bins = 201;

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


