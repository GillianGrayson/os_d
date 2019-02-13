clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'qjx_results';

data_path = sprintf('%s/%s', data_path, prefix);

sys_id = 1;
task_id = 1;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 100;
num_tp_periods = 1337;
num_obs_periods = 10000;
ex_deep = 16;
rk_ns = 10000;

T = 4.8;
jcs_drv_ampl = 4.8;

N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;
jcs_drv_part_1 = 0.98 * T;
jcs_drv_part_2 = 1.0 * T;
jcs_prm_alpha = 5;
start_type = 0;
start_state = 0;

sys_size = N;

num_runs = 1;


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

for run_id = 1:num_runs
    
    ss = (run_id - 1) * num_trajectories;
    
    path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/N_%d/diss_%d_%0.4f_%0.4f/drv_%0.4f_%0.4f_%0.4f/prm_%0.4f/start_%d_%d/ss_%d', ...
        data_path, ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        ex_deep, ...
        rk_ns, ...
        num_tp_periods, ...
        num_obs_periods, ...
        N, ...
        diss_type, ...
        diss_gamma, ...
        diss_phase, ...
        jcs_drv_part_1, ...
        jcs_drv_part_2, ...
        jcs_drv_ampl, ...
        jcs_prm_alpha, ...
        start_type, ...
        start_state, ...
        ss);
    
    suffix = sprintf('config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)', ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        ss, ...
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
    
    for trajectory_id = 1:num_trajectories
        
        path = sprintf('%s/jump_times_%d_%s.txt', path_to_folder, trajectory_id - 1, suffix);
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

jcs_drv_ampl = jcs_drv_ampl

non_inc_count = non_inc_count
norm = sum(curr_pdf)

for bin_id = 1 : num_bins
    curr_pdf(bin_id) = curr_pdf(bin_id) / (norm * bin_diff(bin_id));
end

norm_check = 0;
for bin_id = 1 : num_bins
    norm_check = norm_check + curr_pdf(bin_id) * bin_diff(bin_id);
end
norm_check = norm_check

suffix = sprintf('config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    ss, ...
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
    start_state)
	
[x_min, x_max] = find_range(bin_centers, curr_pdf)
[alpha, coef, R2, yy, R2_, cnt] = powerlaw_regression(bin_centers, curr_pdf, x_min, x_max)
left = find(bin_centers==x_min)
right = find(bin_centers==x_max)

fig = figure;
hLine = plot(bin_centers, curr_pdf);
hold all;
hLine = plot(bin_centers(left - 1:right + 1), yy(left - 1:right + 1));
set(gca, 'FontSize', 30);
xlabel('$f_0$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\Delta t$', 'Interpreter', 'latex');
set(gca,'XScale','log');
set(gca,'YScale','log');

savefig(sprintf('%s/times_between_jumps_%s.fig', home_figures_path, suffix));

close(fig);