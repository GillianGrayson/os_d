clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 3;
task_id = 7;
prop_id = 0;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 1000;
num_obs_periods = 1000;

Nc = 8;

W_seed_s = 1;
W_seed_num = 100;
W_seed_f = W_seed_s + W_seed_num - 1;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

W = 20;
U = 1;
J = 1;

start_type = 0;
start_state = 49;

lpn_type = 0;
lpn_delta_f = 0.001;
lpn_delta_s = 0.000001;

ss = 0;
mns = 1000000;

num_trajectories = 200;
num_target_trajectories = num_trajectories / 2;
num_runs = 1;

obs_min = -1;
ons_max = 100;
obs_num_bins = 400;
obs_shift = (ons_max - obs_min) / obs_num_bins;
obs_bins = linspace(obs_min + 0.5 * obs_shift, ons_max - 0.5 * obs_shift, obs_num_bins);


obs_pdf_deep = zeros(obs_num_bins, 1);

non_inc_deep = 0;
for W_seed = W_seed_s : W_seed_f
    
    for run_id = 1:num_runs
        
        ss = (run_id - 1) * num_trajectories;
        
        path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/Nс_%d/rnd_%d_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/lpn_%d_%0.4f_%0.4f/ss_%d', ...
            data_path, ...
            sys_id, ...
            task_id, ...
            prop_id, ...
            ex_deep, ...
            rk_ns, ...
            num_tp_periods, ...
            num_obs_periods, ...
            Nc, ...
            W_seed, ...
            W_mns, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J, ...
            start_type, ...
            start_state, ...
            lpn_type, ...
            log10(lpn_delta_s), ...
            log10(lpn_delta_f), ...
            ss);
        
        suffix = sprintf('rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f)', ...
            ss, ...
            mns, ...
            Nc, ...
            W_seed, ...
            W_mns, ...
            diss_type, ...
            diss_phase, ...
            diss_gamma, ...
            W, ...
            U, ...
            J, ...
            start_type, ...
            start_state, ...
            lpn_type, ...
            log10(lpn_delta_s), ...
            log10(lpn_delta_f));
        
        fn = sprintf('%s/num_renorms_%s.txt', path_to_folder, suffix);
        obs_data = importdata(fn);
        curr_obs = obs_data(num_target_trajectories + 1 : end);
        
        for l_id = 1: size(curr_obs, 1)
            curr_l = curr_obs(l_id);
            if curr_l < ons_max && curr_l >= obs_min
                curr_id = floor((curr_l - obs_min) / (ons_max - obs_min + eps) * obs_num_bins) + 1;
                obs_pdf_deep(curr_id) = obs_pdf_deep(curr_id) + 1;
            else
                non_inc_deep = non_inc_deep + 1;
            end
        end
    end
end

non_inc_deep = non_inc_deep
sum_lambda_pdf_deep = sum(obs_pdf_deep);
obs_pdf_deep = obs_pdf_deep / (sum_lambda_pdf_deep * obs_shift);
norm = sum(obs_pdf_deep) * obs_shift;
norm_diff_deep = 1.0 - norm

suffix_save =  sprintf('Nc(%d)_W_seed(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f)', ...
    Nc, ...
    W_seed_s, ...
    W_seed_f, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    W, ...
    U, ...
    J, ...
    start_type, ...
    start_state, ...
    lpn_type, ...
    log10(lpn_delta_s), ...
    log10(lpn_delta_f));


fig = figure;
hLine = plot(obs_bins, obs_pdf_deep');
set(gca, 'FontSize', 30);
xlabel('number of renorms', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$PDF$', 'Interpreter', 'latex');

savefig(sprintf('%s/num_renorms_deep_pdf_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/num_renorms_deep_pdf_%s.pdf', home_figures_path, suffix_save));

close(fig)