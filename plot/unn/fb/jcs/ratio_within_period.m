clear all;

path = '../../../../data/jcs';



task = 1;
N = 200;
drv_ampl_start = 0.05;
drv_ampl_shift = 0.05;
drv_ampl_num = 100;

dt = 1;
dp = 0;
g  = 0.1;
g_add = 0.1;

prm_alpha = 5.0;

num_points = 201;

num_first_skips = floor(N/3);
num_last_skips = floor(N/3);

ampls = zeros(drv_ampl_num, 1);
ratio_avg = zeros(drv_ampl_num, 1);
for ampl_id = 1:drv_ampl_num
    ampl = drv_ampl_start + (ampl_id - 1) * drv_ampl_shift
    ampls(ampl_id) = ampl;
    
    fn = sprintf('%s/main_%d/N_%d/prm_%0.4f_%0.4f/diss_%d_%0.4f_%0.4f_%0.4f/evals.txt', ...
        path, ...
        task, ...
        N, ...
        ampl, ...
        prm_alpha, ...
        dt, ...
        dp, ...
        g, ...
        g_add);
    all_evals = importdata(fn);
    
    ratio = 0;
    for p_id = 1:num_points
        evals = all_evals((p_id - 1) * N + 1 : p_id * N);
        evals = sort(evals);
        
        cutted_evals = evals(num_first_skips + 1:size(evals, 1) - num_last_skips);
        size_cutted = size(cutted_evals, 1);
        
        s_n = zeros(size_cutted-1, 1);
        r_n = zeros(size_cutted-2, 1);
        
        for s_id = 1 : (size_cutted-1)
            s_n(s_id) = cutted_evals(s_id + 1) - cutted_evals(s_id);
        end
        
        for s_id = 1 : (size_cutted-2)
            r_n(s_id) = s_n(s_id + 1) / s_n(s_id);
        end
        
        for s_id = 1 : (size_cutted-2)
            r_n(s_id) = min(r_n(s_id), 1.0/r_n(s_id));
        end
        
        r_n_avg = 0.0;
        for s_id = 1 : (size_cutted-2)
            r_n_avg = r_n_avg + r_n(s_id);
        end
        
        r_n_avg = r_n_avg / (size_cutted-2);
        
        ratio = ratio + r_n_avg;
    end
    ratio = ratio / num_points;
    
    ratio_avg(ampl_id) = ratio;

end

fig = figure;
hLine = plot(ampls, ratio_avg, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$f_0$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$ratio$', 'Interpreter', 'latex');
hold all;
propertyeditor(fig)

