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
num_obs_periods = 1000;
ex_deep = 16;
rk_ns = 10000;

N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;
jcs_drv_part_1 = 0.98;
jcs_drv_part_2 = 1;
jcs_drv_ampl = 0.1;
jcs_prm_alpha = 5;
start_type = 0;
start_state = 0;

sys_size = N;

num_runs = 1;

ampl_begin = 0.025;
ampl_step = 0.025;
ampl_num = 20;
ampls = zeros(ampl_num, 1);

bin_begin = 1e-10;
num_decades = 15;
bin_end = bin_begin * 10.^num_decades;
num_bin_per_decade = 5;

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

global_pdf = zeros(ampl_num, num_bins);

for ampl_id = 1:ampl_num
    
    ampl_id = ampl_id
    ampl = ampl_begin + ampl_step * (ampl_id - 1);
    ampls(ampl_id) = ampl;
    
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
            ampl, ...
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
            ampl, ...
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
    
    global_pdf(ampl_id, :) = curr_pdf;
    
end

global_pdf = global_pdf';

suffix = sprintf('config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_var)_prm(%0.4f)_start(%d_%d)', ...
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
    jcs_prm_alpha, ...
    start_type, ...
    start_state);

fig = figure;
hLine = imagesc(ampls, log10(bin_centers), global_pdf);
set(gca, 'FontSize', 30);
xlabel('$f_0$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\Delta t$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF');
set(gca,'YDir','normal');

savefig(sprintf('%s/times_between_jumps_from_ampl_%s.fig', home_figures_path, suffix));

close(fig);