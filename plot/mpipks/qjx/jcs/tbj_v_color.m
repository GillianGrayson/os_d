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
num_trajectories = 10;
num_tp_periods = 1337;
num_obs_periods = 10000;
ex_deep = 16;
rk_ns = 10000;

T = 2.0;
jcs_prm_alpha = 5.0;

N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;
jcs_drv_part_1 = 0.98 * T;
jcs_drv_part_2 = 1 * T;
jcs_drv_ampl = 0.1;
start_type = 0;
start_state = 0;

num_runs = 1;

ampl_begin = 0.05;
ampl_step = 0.05;
ampl_num = 100;
ampls = linspace(ampl_begin, ampl_num * ampl_step, ampl_num);

num_bins = 1000;
v_begin = 0;
v_end = 10;
v_step = (v_end - v_begin) / num_bins;
v_bins = linspace(v_begin +  0.5 * v_step, v_end - 0.5 * v_step, num_bins);

non_inc_count = 0;

pdf = zeros(ampl_num, num_bins);



for ampl_id = 1:ampl_num
    
    ampl = ampl_begin + ampl_step * (ampl_id - 1)
    
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
            jcs_drv_part_1 * T, ...
            jcs_drv_part_2 * T, ...
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
            jcs_drv_part_1 * T, ...
            jcs_drv_part_2 * T, ...
            ampl, ...
            jcs_prm_alpha, ...
            start_type, ...
            start_state);
        
        for trajectory_id = 1:num_trajectories
            
            path = sprintf('%s/jump_times_%d_%s.txt', path_to_folder, trajectory_id - 1, suffix);
            jump_times = importdata(path);
            
            path = sprintf('%s/jump_etas_%d_%s.txt', path_to_folder, trajectory_id - 1, suffix);
            jump_etas = importdata(path);
            
            num_sections = size(jump_times, 1);
            
            for section_id = 2:num_sections
                
                curr_time = jump_times(section_id) - jump_times(section_id - 1);
                curr_dist = -log10(jump_etas(section_id));
                curr_v = curr_dist / curr_time;
                
                if ((curr_v >= v_begin) && (curr_v <= v_end))
                    id = floor((curr_v - v_begin) * num_bins / (v_end - v_begin + 1e-8)) + 1;
                    curr_pdf(id) = curr_pdf(id) + 1;
                else
                    non_inc_count = non_inc_count + 1;
                end
                
            end
        end    
    end
    
    
    curr_pdf = curr_pdf / (sum(curr_pdf) * v_step);
    norm = sum(curr_pdf) * v_step
    non_inc_count = non_inc_count  
    
    pdf(ampl_id, :) = curr_pdf / max(curr_pdf);
end



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
hLine = imagesc(ampls, v_bins, pdf');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$V$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$PDF$', 'Interpreter', 'latex');
set(gca,'YDir','normal');
savefig(sprintf('%s/v_color_%s.fig', home_figures_path, suffix));

