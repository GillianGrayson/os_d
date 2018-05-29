clear all;

path = '../../../../data/jcs';

task = 2;
N = 10;
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
max_mult = zeros(drv_ampl_num, 1);
for ampl_id = 1:drv_ampl_num
    ampl = drv_ampl_start + (ampl_id - 1) * drv_ampl_shift
    ampls(ampl_id) = ampl;
    
    fn = sprintf('%s/main_%d/N_%d/prm_%0.4f_%0.4f/diss_%d_%0.4f_%0.4f_%0.4f/floquet_evals.txt', ...
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
    
    abs_mults = zeros(size(all_evals, 1));
    
    for mult_id = 1:size(all_evals, 1)
        abs_mults(mult_id) = sqrt(all_evals(mult_id, 1) * all_evals(mult_id, 1) + all_evals(mult_id, 2) * all_evals(mult_id, 2));
    end
    
    
    max_mult(ampl_id) = max(abs_mults);

end

fig = figure;
hLine = plot(ampls, max_mult, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$f_0$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$|mult|_{max}$', 'Interpreter', 'latex');
hold all;
propertyeditor(fig)

