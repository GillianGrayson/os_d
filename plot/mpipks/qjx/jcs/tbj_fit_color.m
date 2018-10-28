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

num_runs = 1;

ampl_begin = 0.05;
ampl_step = 0.05;
ampl_num = 100;
ampls = linspace(ampl_begin, ampl_num * ampl_step, ampl_num);

T_begin = 0.05;
T_step = 0.05;
T_num = 100;
Ts = linspace(T_begin, T_num * T_step, T_num);

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

alphas = zeros(ampl_num, T_num);
decades = zeros(ampl_num, T_num);


for T_id = 1:T_num
    
    T = T_begin + T_step * (T_id - 1);
    
    for ampl_id = 1:ampl_num
        
		T = T
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
				
				data = [];
            end
            
        end
        non_inc_count = non_inc_count;
        norm = sum(curr_pdf);
        for bin_id = 1 : num_bins
            curr_pdf(bin_id) = curr_pdf(bin_id) / (norm * bin_diff(bin_id));
        end
        norm_check = 0;
        for bin_id = 1 : num_bins
            norm_check = norm_check + curr_pdf(bin_id) * bin_diff(bin_id);
        end
        norm_check = norm_check;
        
        [x_min, x_max] = find_range(bin_centers, curr_pdf);
        [alpha, coef, R2, yy, R2_, cnt] = powerlaw_regression(bin_centers, curr_pdf, x_min, x_max);
		
		curr_pdf = [];
		
		x_min = x_min;
		x_max = x_max;
		
		alpha = alpha;
		decade = log10(x_max) - log10(x_min);
        
        alphas(ampl_id, T_id) = abs(alpha);
        decades(ampl_id, T_id) = log10(x_max) - log10(x_min);
        
    end
	
	name = sprintf('T_%0.4f.txt', T);
	file_id = fopen(name, 'w');
	for id = 1:ampl_num
		fprintf(file_id, '%0.16e %0.16e\n', alphas(id, T_id), decades(id, T_id));
	end
	fclose(file_id);
	
    
end

suffix = sprintf('config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    ss, ...
    mns, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    jcs_prm_alpha, ...
    start_type, ...
    start_state);


fig = figure;
hLine = imagesc(ampls, Ts, alphas');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$T$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\alpha');
set(gca,'YDir','normal');
savefig(sprintf('%s/alphas_color_%s.fig', home_figures_path, suffix));

fig = figure;
hLine = imagesc(ampls, Ts, decades');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$T$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'decades');
set(gca,'YDir','normal');
savefig(sprintf('%s/decades_color_%s.fig', home_figures_path, suffix));




function [alpha, coef, R2, yy, R2_, cnt] = powerlaw_regression(x, y, mn, mx)
    if nargin < 3
        mn = min(x);
    end
    if nargin < 4
        mx = max(x);
    end
    xx = x;
    y0 = y;
    ids = (y > 0) & (mn < x) & (x < mx);
    x = x(ids);
    y = y(ids);
    
    cnt = size(x);
    
    logx = log(x);
    logy = log(y);
    X = [ones(length(logx),1) logx];
    c = X \ logy;
    logyy = X * c;
    S1 = sum((logy - logyy).^2);
    S2 = sum((logy - mean(logy)).^2);
    R2 = 1 - S1 / S2;
    b = c(1);
    a = c(2);
    coef = exp(b);
    alpha = a;
    yy = coef * xx .^ alpha;
    R2_ = 1 - sum((y0 - yy).^2) / sum((y0 - mean(y0)).^2);
	
end

function [i_max, j_max] = find_range(x, y)
    i_max = 1; j_max = 1;
    ids = y > 0;
    x = x(ids);
    y = y(ids);
    len = 1.0;
    R2_max = 0;
    for i = 1:size(x)
        for j = i:size(x)
            [~, ~, R2_pow, ~, ~, cnt] = powerlaw_regression(x, y, x(i), x(j));
            if cnt < 10
                continue;
            end
            curlen = log10(x(j)) - log10(x(i));
            if R2_pow > 0.98 &&  curlen > 1.0
                 curlen = curlen + 120 * (R2_pow - 0.98);
                if abs(len - curlen) <= 0.01
                    if R2_pow > R2_max
                        len = curlen;
                        R2_max = R2_pow;
                        i_max = x(i);
                        j_max = x(j);
                    end
                elseif len + 1 < curlen
                    len = curlen;
                    R2_max = R2_pow;
                    i_max = x(i);
                    j_max = x(j);
                end
            end
        end
    end 
end

