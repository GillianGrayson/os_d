clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';

data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';


sys_id = 2; 
task_id = 1; 
prop_id = 0; 
seed = 0; 
mns = 1000000;
num_trajectories = 10;
num_tp_periods = 100;
num_obs_periods = 100; 
ex_deep = 16;
rk_ns = 10000;

ps_drv_ampl = 3.25;
T = 2.0;
g = 1.0;
d_start = 0.0;
d_shift = 0.1;

diss_type = 1;
ps_diss_w = 0.05;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
ps_drv_part_1 = 1.00 * T;
ps_drv_part_2 = 1.00 * T;
ps_prm_alpha = 5;
start_type = 0;
start_state = 0;

num_runs = 10;

num_bins = 200;

num_points = 100;
params = zeros(num_points, 1);
theta_re = zeros(num_points, num_bins);
theta_im = zeros(num_points, num_bins);
J_plus_re = zeros(num_points, num_bins);
J_plus_im = zeros(num_points, num_bins);
J_z = zeros(num_points, num_bins);

min_theta_re = -3;
max_theta_re = 3;
shift_theta_re = (max_theta_re - min_theta_re) / num_bins;
int_theta_re = zeros(num_bins, 1);
for int_id = 1:num_bins
    int_theta_re(int_id) = min_theta_re + shift_theta_re * int_id - 0.5 * shift_theta_re;
end

min_theta_im = -3;
max_theta_im = 3;
shift_theta_im = (max_theta_im - min_theta_im) / num_bins;
int_theta_im = zeros(num_bins, 1);
for int_id = 1:num_bins
    int_theta_im(int_id) = min_theta_im + shift_theta_im * int_id - 0.5 * shift_theta_im;
end

min_J_plus_re = -ps_num_spins * 0.5 - 0.1;
max_J_plus_re = ps_num_spins * 0.5 + 0.1;
shift_J_plus_re = (max_J_plus_re - min_J_plus_re) / num_bins;
int_J_plus_re = zeros(num_bins, 1);
for int_id = 1:num_bins
    int_J_plus_re(int_id) = min_J_plus_re + shift_J_plus_re * int_id - 0.5 * shift_J_plus_re;
end

min_J_plus_im = -ps_num_spins * 0.5 - 0.1;
max_J_plus_im = ps_num_spins * 0.5 + 0.1;
shift_J_plus_im = (max_J_plus_im - min_J_plus_im) / num_bins;
int_J_plus_im = zeros(num_bins, 1);
for int_id = 1:num_bins
    int_J_plus_im(int_id) = min_J_plus_im + shift_J_plus_im * int_id - 0.5 * shift_J_plus_im;
end

min_J_z = -ps_num_spins * 0.5 - 0.1;
max_J_z = ps_num_spins * 0.5 + 0.1;
shift_J_z = (max_J_z - min_J_z) / num_bins;
int_J_z = zeros(num_bins, 1);
for int_id = 1:num_bins
    int_J_z(int_id) = min_J_z + shift_J_z * int_id - 0.5 * shift_J_z;
end



for param_id = 1:num_points
    
    ps_prm_d = d_start + (param_id - 1) * d_shift
    params(param_id) = ps_prm_d;
	ps_prm_g = g;
    
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
        data_spec = importdata(path);
        
        path = sprintf('%s/spec_2_evo_%s.txt', path_to_folder, suffix);
        data_spec_2 = importdata(path);
        
        path = sprintf('%s/spec_3_evo_%s.txt', path_to_folder, suffix);
        data_spec_3 = importdata(path);
        
        theta_re_curr = data_spec(:, 1);
        theta_im_curr = data_spec(:, 2);
        J_plus_re_curr = data_spec_3(:, 1);
        J_plus_im_curr = data_spec_3(:, 2);
        J_z_curr = data_spec_2(:, 1);
        
        for tr_id = 1:num_trajectories-1
            theta_re_curr = vertcat(theta_re_curr, data_spec(:, 1 + 2*tr_id));
            theta_im_curr = vertcat(theta_im_curr, data_spec(:, 2 + 2*tr_id));
            J_plus_re_curr = vertcat(J_plus_re_curr, data_spec_3(:, 1 + 2*tr_id));
            J_plus_im_curr = vertcat(J_plus_im_curr, data_spec_3(:, 2 + 2*tr_id));
            J_z_curr = vertcat(J_z_curr, data_spec_2(:, 1 + 2*tr_id));
        end
        
        global_size = size(J_z_curr, 1);
        
        for d_id = 1:global_size

            int_theta_re_id = floor((theta_re_curr(d_id) - min_theta_re) * num_bins / (max_theta_re - min_theta_re + 10e-8)) + 1;
            int_theta_im_id = floor((theta_im_curr(d_id) - min_theta_im) * num_bins / (max_theta_im - min_theta_im + 10e-8)) + 1;
            int_J_plus_re_id = floor((J_plus_re_curr(d_id) - min_J_plus_re) * num_bins / (max_J_plus_re - min_J_plus_re + 10e-8)) + 1;
            int_J_plus_im_id = floor((J_plus_im_curr(d_id) - min_J_plus_im) * num_bins / (max_J_plus_im - min_J_plus_im + 10e-8)) + 1;
            int_J_z_id = floor((J_z_curr(d_id) - min_J_z) * num_bins / (max_J_z - min_J_z + 10e-8)) + 1;
            
            theta_re(param_id, int_theta_re_id) = theta_re(param_id, int_theta_re_id) + 1;
            theta_im(param_id, int_theta_im_id) = theta_im(param_id, int_theta_im_id) + 1;
            J_plus_re(param_id, int_J_plus_re_id) = J_plus_re(param_id, int_J_plus_re_id) + 1;
            J_plus_im(param_id, int_J_plus_im_id) = J_plus_im(param_id, int_J_plus_im_id) + 1;
            J_z(param_id, int_J_z_id) = J_z(param_id, int_J_z_id) + 1;
        end
        
    end
    
    theta_re(param_id, :) = theta_re(param_id, :) / max(theta_re(param_id, :));
    theta_im(param_id, :) = theta_im(param_id, :) / max(theta_im(param_id, :));
    J_plus_re(param_id, :) = J_plus_re(param_id, :) / max(J_plus_re(param_id, :));
    J_plus_im(param_id, :) = J_plus_im(param_id, :) / max(J_plus_im(param_id, :));
    J_z(param_id, :) = J_z(param_id, :) / max(J_z(param_id, :));
    
end

fig = figure;

subplot(5,1,1);
hLine = imagesc(params, int_theta_re, theta_re');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel('$Re(\theta)$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 0.78 0.8 0.16]);
set(h, 'Position', [0.91 0.78 0.02 0.16]);
set(h,'ytick', [0.5 1])
hold all;

subplot(5,1,2);
hLine = imagesc(params, int_theta_im, theta_im');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel('$Im(\theta)$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 0.61 0.8 0.16]);
set(h, 'Position', [0.91 0.61 0.02 0.16]);
set(h,'ytick', [0.5 1])
hold all;

subplot(5,1,3);
hLine = imagesc(params, int_J_plus_re, J_plus_re');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel('$Re(J_{+})$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 0.44 0.8 0.16]);
set(h, 'Position', [0.91 0.44 0.02 0.16]);
set(h,'ytick', [0.5 1])
hold all;

subplot(5,1,4);
hLine = imagesc(params, int_J_plus_im, J_plus_im');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel('$Im(J_{+})$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 0.27 0.8 0.16]);
set(h, 'Position', [0.91 0.27 0.02 0.16]);
set(h,'ytick', [0.5 1])
hold all;

subplot(5,1,5);
hLine = imagesc(params, int_J_z, J_z');
set(gca, 'FontSize', 20);
xlabel('$\delta$', 'Interpreter', 'latex');
ylabel('$J_{z}$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 0.1 0.8 0.16]);
set(h, 'Position', [0.91 0.1 0.02 0.16]);
set(h,'ytick', [0.5 1])
hold all;

suffix_save = sprintf("s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_var_%0.4f)", ...
	ps_num_spins, ...
	ps_num_photons_states, ...
	diss_type, ...
	ps_diss_w, ...
	ps_drv_part_1, ...
	ps_drv_part_2, ...
	ps_drv_ampl, ...
	ps_prm_alpha, ...
    g);

savefig(sprintf('%s/observables_from_d_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/observables_from_d_%s.pdf', home_figures_path, suffix_save));
