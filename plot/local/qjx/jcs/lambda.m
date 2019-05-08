clear all;

sys_id = 1;
task_id = 1;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 10;

N = 200;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
jcs_drv_part_1 = 0.98;
jcs_drv_part_2 = 1.0;
jcs_drv_ampl = 3.2;
jcs_prm_alpha = 5.0;
start_type = 0;
start_state = 0;

path = "../../../../source/cpp/QJX/QJX";

num_points = 100;

lambdas = zeros(num_points, 1);
params = zeros(num_points, 1);

for param_id = 1:num_points
    
    jcs_drv_ampl = 0.05 + (param_id - 1) * 0.05;
    params(param_id) = jcs_drv_ampl;

    suffix = sprintf("rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)", ...
        seed, ...
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

    fn = sprintf('%s/lambda_%s.txt', path, suffix);
    data = importdata(fn);
    
    lambdas(param_id) = mean(data(2:end));

end

fig = figure;
hLine = plot(params, lambdas);
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\lambda$', 'Interpreter', 'latex');
hold all;

propertyeditor('on')


