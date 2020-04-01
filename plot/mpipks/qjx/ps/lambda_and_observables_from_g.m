clear all;

addpath('/home/ivanchen/yusipov/os_lnd/source/matlab/lib')

home_figures_path = '/home/ivanchen/yusipov/os_d/figures';

data_path = '/data/condmat/ivanchen/yusipov/os_d/qjx';

task_id_lambda = 7;
num_trajectories_lambda = 200;
num_runs_lambda = 1;

sys_id = 2;
task_id = 1;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 100;
num_tp_periods = 10;
num_obs_periods = 10;
ex_deep = 16;
rk_ns = 10000;

num_random_obs = 1;
random_obs_seed = 100;
random_obs_type = 2;

lpn_type = -1;
lpn_delta_s = log10(1.0e-4);
lpn_delta_f_high = log10(1.0e-4);
lpn_delta_f_low = log10(1.0e-4);

ampl = 2.75;
T = 4.0;
d = 1.0;
g_start = 0.01;
g_shift = 0.01;
g_num = 1000;
gs = linspace(g_start, g_start + (g_num - 1) * g_shift, g_num);

diss_type = 1;
ps_diss_w = 0.05;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 300;
ps_drv_part_1 = 1.00 * T;
ps_drv_part_2 = 1.00 * T;
ps_prm_alpha = 5;
start_type = 0;
start_state = 0;

num_runs = 1;

num_bins = 200;

lambdas = zeros(g_num, 1);

theta_re_pdf.xs = gs;
theta_re_pdf.x_num_points = g_num;
theta_re_pdf.y_num_bins = num_bins;
theta_re_pdf.x_label = '$g$';
theta_re_pdf.y_label = '$Re(\theta)$';

theta_im_pdf.xs = gs;
theta_im_pdf.x_num_points = g_num;
theta_im_pdf.y_num_bins = num_bins;
theta_im_pdf.x_label = '$g$';
theta_im_pdf.y_label = '$Im(\theta)$';

J_plus_re_pdf.xs = gs;
J_plus_re_pdf.x_num_points = g_num;
J_plus_re_pdf.y_num_bins = num_bins;
J_plus_re_pdf.x_label = '$g$';
J_plus_re_pdf.y_label = '$Re(\nu)$';

J_plus_im_pdf.xs = gs;
J_plus_im_pdf.x_num_points = g_num;
J_plus_im_pdf.y_num_bins = num_bins;
J_plus_im_pdf.x_label = '$g$';
J_plus_im_pdf.y_label = '$Im(\nu)$';

J_z_pdf.xs = gs;
J_z_pdf.x_num_points = g_num;
J_z_pdf.y_num_bins = num_bins;
J_z_pdf.x_label = '$g$';
J_z_pdf.y_label = '$\eta$';

num_global = (num_obs_periods + 1) * num_trajectories * num_runs;

theta_re = zeros(g_num, num_global);
theta_im = zeros(g_num, num_global);
J_plus_re = zeros(g_num, num_global);
J_plus_im = zeros(g_num, num_global);
J_z = zeros(g_num, num_global);

for g_id = 1:g_num
    
    g = gs(g_id);
	fprintf('g = %0.16e\n', g);
    
    for run_id = 1:num_runs_lambda
        
        ss = (run_id - 1) * num_trajectories_lambda;
        
        path_to_folder = sprintf('%s/main_%d_%d_%d/lpn_%d_%0.4f_%0.4f_%0.4f/run_%d_%d_%d_%d/obs_%d_%d_%d/N_%d_%d/diss_%d_%0.4f/drv_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
            data_path, ...
            sys_id, ...
            task_id_lambda, ...
            prop_id, ...
            lpn_type, ...
            lpn_delta_s, ...
            lpn_delta_f_high, ...
            lpn_delta_f_low, ...
            ex_deep, ...
            rk_ns, ...
            num_tp_periods, ...
            num_obs_periods, ...
            num_random_obs, ...
            random_obs_seed, ...
            random_obs_type, ...
            ps_num_spins, ...
            ps_num_photons_states, ...
            diss_type, ...
            ps_diss_w, ...
            ps_drv_part_1, ...
            ps_drv_part_2, ...
            ampl, ...
            ps_prm_alpha, ...
            d, ...
            g, ...
            start_type, ...
            start_state, ...
            ss);
        
        suffix = sprintf("setup(%d_%d_%d)_rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)", ...
            sys_id, ...
            task_id_lambda, ...
            prop_id, ...
            ss, ...
            mns, ...
            ps_num_spins, ...
            ps_num_photons_states, ...
            diss_type, ...
            ps_diss_w, ...
            ps_drv_part_1, ...
            ps_drv_part_2, ...
            ampl, ...
            ps_prm_alpha, ...
            d, ...
            g, ...
            start_type, ...
            start_state, ...
            lpn_type, ...
            lpn_delta_s, ...
            lpn_delta_f_high, ...
            lpn_delta_f_low);
        
        path = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
        data = importdata(path);
        
        lambdas(g_id) = lambdas(g_id) + mean(data(num_trajectories_lambda / 2 + 1:end));
    end
    
    for run_id = 1:num_runs
        
        ss = (run_id - 1) * num_trajectories;
        
        path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/obs_%d_%d_%d/N_%d_%d/diss_%d_%0.4f/drv_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
            data_path, ...
            sys_id, ...
            task_id, ...
            prop_id, ...
            ex_deep, ...
            rk_ns, ...
            num_tp_periods, ...
            num_obs_periods, ...
            num_random_obs, ...
            random_obs_seed, ...
            random_obs_type, ...
            ps_num_spins, ...
            ps_num_photons_states, ...
            diss_type, ...
            ps_diss_w, ...
            ps_drv_part_1, ...
            ps_drv_part_2, ...
            ampl, ...
            ps_prm_alpha, ...
            d, ...
            g, ...
            start_type, ...
            start_state, ...
            ss);
        
        suffix = sprintf("setup(%d_%d_%d)_rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
            sys_id, ...
            task_id, ...
            prop_id, ...
            ss, ...
            mns, ...
            ps_num_spins, ...
            ps_num_photons_states, ...
            diss_type, ...
            ps_diss_w, ...
            ps_drv_part_1, ...
            ps_drv_part_2, ...
            ampl, ...
            ps_prm_alpha, ...
            d, ...
            g, ...
            start_type, ...
            start_state);
        
        path = sprintf('%s/spec_evo_%s.txt', path_to_folder, suffix);
        data_spec = importdata(path);
        
        path = sprintf('%s/spec_2_evo_%s.txt', path_to_folder, suffix);
        data_spec_2 = importdata(path);
        
        path = sprintf('%s/spec_3_evo_%s.txt', path_to_folder, suffix);
        data_spec_3 = importdata(path);
        
        theta_re_curr = data_spec(:, 1:2:end);
        theta_re_curr = theta_re_curr(:);
        theta_im_curr = data_spec(:, 2:2:end);
        theta_im_curr = theta_im_curr(:);
        J_plus_re_curr = data_spec_3(:, 1:2:end);
        J_plus_re_curr = J_plus_re_curr(:);
        J_plus_im_curr = data_spec_3(:, 2:2:end);
        J_plus_im_curr = J_plus_im_curr(:);
        J_z_curr = data_spec_2(:, 1:2:end);
        J_z_curr = J_z_curr(:);
        
        s_id = (run_id - 1) * (num_obs_periods + 1) * num_trajectories + 1;
        f_id = run_id * (num_obs_periods + 1) * num_trajectories;
        
        theta_re(g_id, s_id:f_id) = theta_re_curr;
        theta_im(g_id, s_id:f_id) = theta_im_curr;
        J_plus_re(g_id, s_id:f_id) = J_plus_re_curr;
        J_plus_im(g_id, s_id:f_id) = J_plus_im_curr;
        J_z(g_id, s_id:f_id) = J_z_curr;
        
    end
end


theta_re_pdf.y_bin_s = min(theta_re, [], 'all');
theta_re_pdf.y_bin_f = max(theta_re, [], 'all');
theta_re_pdf = oqs_pdf_2d_lead_by_x_setup(theta_re_pdf);
theta_re_pdf = oqs_pdf_2d_lead_by_x_update(theta_re_pdf, theta_re);
theta_re_pdf = oqs_pdf_2d_lead_by_x_release(theta_re_pdf);

theta_im_pdf.y_bin_s = min(theta_im, [], 'all');
theta_im_pdf.y_bin_f = max(theta_im, [], 'all');
theta_im_pdf = oqs_pdf_2d_lead_by_x_setup(theta_im_pdf);
theta_im_pdf = oqs_pdf_2d_lead_by_x_update(theta_im_pdf, theta_im);
theta_im_pdf = oqs_pdf_2d_lead_by_x_release(theta_im_pdf);

J_plus_re_pdf.y_bin_s = min(J_plus_re, [], 'all');
J_plus_re_pdf.y_bin_f = max(J_plus_re, [], 'all');
J_plus_re_pdf = oqs_pdf_2d_lead_by_x_setup(J_plus_re_pdf);
J_plus_re_pdf = oqs_pdf_2d_lead_by_x_update(J_plus_re_pdf, J_plus_re);
J_plus_re_pdf = oqs_pdf_2d_lead_by_x_release(J_plus_re_pdf);

J_plus_im_pdf.y_bin_s = min(J_plus_im, [], 'all');
J_plus_im_pdf.y_bin_f = max(J_plus_im, [], 'all');
J_plus_im_pdf = oqs_pdf_2d_lead_by_x_setup(J_plus_im_pdf);
J_plus_im_pdf = oqs_pdf_2d_lead_by_x_update(J_plus_im_pdf, J_plus_im);
J_plus_im_pdf = oqs_pdf_2d_lead_by_x_release(J_plus_im_pdf);

J_z_pdf.y_bin_s = min(J_z, [], 'all');
J_z_pdf.y_bin_f = max(J_z, [], 'all');
J_z_pdf = oqs_pdf_2d_lead_by_x_setup(J_z_pdf);
J_z_pdf = oqs_pdf_2d_lead_by_x_update(J_z_pdf, J_z);
J_z_pdf = oqs_pdf_2d_lead_by_x_release(J_z_pdf);

fig = figure;

shift_y = 0.85/6;
starts_y = zeros(6, 1);
for i = 1:6
    starts_y(i) = 0.1 + (i-1) * shift_y;
end

subplot(6,1,1);
hLine = imagesc(theta_re_pdf.xs, theta_re_pdf.y_bin_centers, theta_re_pdf.pdf');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel(theta_re_pdf.y_label, 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 starts_y(6) 0.8 0.14]);
set(h, 'Position', [0.91 starts_y(2) 0.02 0.7085]);
set(h,'ytick', [0.25 0.5 0.75 1])
hold all;

subplot(6,1,2);
hLine = imagesc(theta_im_pdf.xs, theta_im_pdf.y_bin_centers, theta_im_pdf.pdf');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel(theta_im_pdf.y_label, 'Interpreter', 'latex');
colorbar('off')
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 starts_y(5) 0.8 0.14]);
hold all;

subplot(6,1,3);
hLine = imagesc(J_plus_re_pdf.xs, J_plus_re_pdf.y_bin_centers, J_plus_re_pdf.pdf');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel(J_plus_re_pdf.y_label, 'Interpreter', 'latex');
colorbar('off')
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 starts_y(4) 0.8 0.14]);
hold all;

subplot(6,1,4);
hLine = imagesc(J_plus_im_pdf.xs, J_plus_im_pdf.y_bin_centers, J_plus_im_pdf.pdf');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel(J_plus_im_pdf.y_label, 'Interpreter', 'latex');
colorbar('off')
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 starts_y(3) 0.8 0.14]);
hold all;

subplot(6,1,5);
hLine = imagesc(J_z_pdf.xs, J_z_pdf.y_bin_centers, J_z_pdf.pdf');
set(gca, 'FontSize', 20);
xlabel('')
set(gca,'xticklabel',{[]})
ylabel(J_z_pdf.y_label, 'Interpreter', 'latex');
colorbar('off')
set(gca, 'FontSize', 20);
set(gca,'YDir','normal');
set(gca, 'Position', [0.1 starts_y(2) 0.8 0.14]);
hold all;

subplot(6,1,6);
hLine = plot(gs, lambdas, 'LineWidth', 2);
hold all;
hLine = plot([gs(1) gs(end)], [0 0], 'k', 'LineStyle','--');
set(gca, 'FontSize', 20);
xlabel('$g$', 'Interpreter', 'latex');
ylabel('$\lambda$', 'Interpreter', 'latex');
set(gca, 'Position', [0.1 starts_y(1) 0.8 0.14]);
hold all;

suffix_save = sprintf("s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_var)", ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    diss_type, ...
    ps_diss_w, ...
    ps_drv_part_1, ...
    ps_drv_part_2, ...
    ampl, ...
    ps_prm_alpha, ...
    d);

savefig(sprintf('%s/lambda_and_observables_from_g_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_and_observables_from_g_%s.pdf', home_figures_path, suffix_save));

