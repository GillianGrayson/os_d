clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 3;
task_id = 1;
prop_id = 0;

T = 10;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 500;
num_obs_periods = 100000;

Nc = 8;

obs = 'obs_1_100_2';

W_seed_s = 1;
W_seed_num = 2;
W_seed_f = W_seed_s + W_seed_num - 1;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

U = 1;
J = 1;

start_type = 0;
start_state = 0;

ss = 0;
mns = 1000000;

num_trajectories = 50;
num_runs = 1;

W_begin = 1.0;
W_step = 0.2;
W_num = 1;
Ws = zeros(W_num, 1);

for W_id = 1:W_num
    
    W = W_begin + W_step * (W_id - 1);
    fprintf('W = %0.16e\n', W);
    Ws(W_id) = W;
    
    bin_begin = 1e-10;
    num_decades = 15;
    bin_end = bin_begin * 10.^num_decades;
    num_bin_per_decade = 10;
    num_bins = num_bin_per_decade * num_decades;
    bin_borders = zeros(num_bins + 1, 1);
    bin_centers = zeros(num_bins, 1);
    for bin_id = 1 : num_bins + 1
        bin_borders(bin_id) = bin_begin * 10.^((bin_id - 1) / num_bin_per_decade);
        if (bin_id <= num_bins)
            bin_centers(bin_id) = bin_begin * 10.^((bin_id - 1 + 0.5) / num_bin_per_decade);
        end
    end
    bin_diff = diff(bin_borders);
    non_inc_count = 0;
    curr_pdf = zeros(num_bins, 1);
    
    
    for W_seed = W_seed_s : W_seed_f
        
        for run_id = 1:num_runs
            
            ss = (run_id - 1) * num_trajectories;
        
            path_to_folder = sprintf('%s/main_%d_%d_%d/T_%0.4f/run_%d_%d_%d_%d/Nс_%d/%s/rnd_%d_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
                data_path, ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                T, ...
                ex_deep, ...
                rk_ns, ...
                num_tp_periods, ...
                num_obs_periods, ...
                Nc, ...
                obs, ...
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
                ss);

            suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
                sys_id, ...
                task_id, ...
                prop_id, ...
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
                start_state);
        
            for trajectory_id = 1:num_trajectories

                path = sprintf('%s/jump_times_%d_%s.txt', path_to_folder, trajectory_id - 1, suffix)
                data = importdata(path);

                curr_diff = diff(data);
                curr_size = size(curr_diff, 1);

                for jump_diff_id = 1:curr_size
                    curr_tmp = curr_diff(jump_diff_id);
                    if ((curr_tmp >= bin_begin) && (curr_tmp <= bin_end))
                        bin_id = floor((log10(curr_tmp) - log10(bin_begin)) * num_bins / (log10(bin_end) - log10(bin_begin) + eps)) + 1;
                        curr_pdf(bin_id) = curr_pdf(bin_id) + 1;
                    else
                        non_inc_count = non_inc_count + 1;
                    end
                end
            end
        end
    end
    
    non_inc_count = non_inc_count
    norm = sum(curr_pdf);
    for bin_id = 1 : num_bins
        curr_pdf(bin_id) = curr_pdf(bin_id) / (norm * bin_diff(bin_id));
    end
    norm_check = 0;
    for bin_id = 1 : num_bins
        norm_check = norm_check + curr_pdf(bin_id) * bin_diff(bin_id);
    end
    norm_check = norm_check
    
    [x_min, x_max] = oqs_find_range(bin_centers, curr_pdf);
    [alpha, coef, R2, yy, R2_, cnt] = oqs_powerlaw_regression(bin_centers, curr_pdf, x_min, x_max);
    left = find(bin_centers==x_min);
    right = find(bin_centers==x_max);
    
    if ~isnan(R2)
        alpha = abs(alpha);
        decades = log10(x_max) - log10(x_min);
    else
        alpha = 0;
        decades = 0;
    end
    
	alpha = alpha
    decades = decades
	
    suffix_save =  sprintf('setup(%d_%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(var_%0.4f_%0.4f)_start(%d_%d)', ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        Nc, ...
        W_seed, ...
        W_mns, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        U, ...
        J, ...
        start_type, ...
        start_state);
    
    fig = figure;
    hLine = plot(bin_centers, curr_pdf);
    hold all;
    hLine = plot(bin_centers(left - 1:right + 1), yy(left - 1:right + 1));
    set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(gca, 'FontSize', 30);
    xlabel('$\Delta t$', 'Interpreter', 'latex');
    set(gca, 'FontSize', 30);
    ylabel('$PDF$', 'Interpreter', 'latex');
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    savefig(sprintf('%s/tbj_%s.fig', home_figures_path, suffix_save));
    h=gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', sprintf('%s/tbj_%s.pdf', home_figures_path, suffix_save));
    
end
