clear all;

seed = 0;
mns = 1000000;
num_trajectories = 10;

diss_type = 1;
ps_num_spins = 1;
ps_num_photons_states = 200;
ps_drv_part_1 = 0.98; 
ps_drv_part_2 = 1.00; 
ps_drv_ampl = 3.2;
ps_prm_alpha = 5;
ps_prm_d = 0.00;
ps_prm_g = 0.00;
ps_diss_w = 0.00;
start_type = 0;
start_state = 0;

ps_num_spins_states = 2^ps_num_spins;



path = "../../../../source/cpp/QJX/QJX";

num_points = 100;

lambdas = zeros(num_points, 1);
params = zeros(num_points, 1);

for param_id = 1:num_points
    
    ps_drv_ampl = 0.05 + (param_id - 1) * 0.05;
    params(param_id) = ps_drv_ampl;

    suffix = sprintf("rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
        seed, ...
        mns, ...
        ps_num_spins, ...
        ps_num_photons_states, ...
        diss_type, ...
        ps_diss_w, ...
        ps_drv_part_1, ...
        ps_drv_part_2, ...
        ps_drv_ampl, ...
        ps_prm_alpha, ...
        ps_prm_d, ...
        ps_prm_g, ...
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