clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';

data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 2;
task_id = 7;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 20;
num_tp_periods = 100;
num_obs_periods = 100;
ex_deep = 16;
rk_ns = 10000;

is_trans = 1;
T = 2;
g = 1;

diss_type = 1;
ps_diss_w = 0.05;
ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
ps_num_photons_states = 200;
ps_drv_part_1 = 1.00 * T;
ps_drv_part_2 = 1.00 * T;
ps_prm_alpha = 5;
start_type = 0;
start_state = 0;

num_runs = 1;

ampl_begin = 0.05;
ampl_step = 0.05;
ampl_num = 100;
ampls = linspace(ampl_begin, ampl_begin + (ampl_num - 1) * ampl_step, ampl_num);

d_begin = 0.0;
d_step = 0.1;
d_num = 100;
ds = linspace(d_begin, d_begin + (d_num - 1) * d_step, d_num);

lambdas = zeros(ampl_num, d_num);

for d_id = 1:d_num
    
    d = d_begin + (d_id - 1) * d_step
    
    for ampl_id = 1:ampl_num
        
        ampl = ampl_begin + (ampl_id - 1) * ampl_step;
        
        for run_id = 1:num_runs
            
            ss = (run_id - 1) * num_trajectories;
            
            path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/N_%d_%d/diss_%d_%0.4f/drv_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
                data_path, ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                ex_deep, ...
                rk_ns, ...
                num_tp_periods, ...
                num_obs_periods, ...
                ps_num_spins, ...
                ps_num_photons_states, ...
                diss_type, ...
                ps_diss_w, ...
                ps_drv_part_1, ...
                ps_drv_part_2, ...
                ampl, ...
                ps_prm_alpha, ...
                d, ...
                g, ...
                start_type, ...
                start_state, ...
                ss);
            
            suffix = sprintf("rnd(%d_%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
                ss, ...
                mns, ...
                ps_num_spins, ...
                ps_num_photons_states, ...
                diss_type, ...
                ps_diss_w, ...
                ps_drv_part_1, ...
                ps_drv_part_2, ...
                ampl, ...
                ps_prm_alpha, ...
                d, ...
                g, ...
                start_type, ...
                start_state);
            
            path = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
            data = importdata(path);
            
            lambdas(ampl_id, d_id) = lambdas(ampl_id, d_id) + mean(data(num_trajectories / 2 + 1:end));
        end
    end  
end

lambdas = lambdas / num_runs;

if is_trans == 1
	fig = figure;
	hLine = imagesc(ds, ampls, lambdas);
	set(gca, 'FontSize', 30);
	xlabel('$d$', 'Interpreter', 'latex');
	set(gca, 'FontSize', 30);
	ylabel('$A$', 'Interpreter', 'latex');
	colormap hot;
	h = colorbar;
	set(gca, 'FontSize', 30);
	title(h, '\lambda');
	set(gca,'YDir','normal');
else
	fig = figure;
	hLine = imagesc(ampls, ds, lambdas');
	set(gca, 'FontSize', 30);
	xlabel('$A$', 'Interpreter', 'latex');
	set(gca, 'FontSize', 30);
	ylabel('$d$', 'Interpreter', 'latex');
	colormap hot;
	h = colorbar;
	set(gca, 'FontSize', 30);
	title(h, '\lambda');
	set(gca,'YDir','normal');
end

fig = figure;
hLine = imagesc(ampls, ds, lambdas');
set(gca, 'FontSize', 30);
xlabel('$A$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$d$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\lambda');
set(gca,'YDir','normal');

suffix = sprintf("trans(%d)_s(%d)_nps(%d)_diss(%d_%0.4f)_drv(%0.4f_%0.4f_var)_prm(%0.4f_var_%0.4f)", ...
	is_trans, ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    diss_type, ...
    ps_diss_w, ...
    ps_drv_part_1, ...
    ps_drv_part_2, ...
    ps_prm_alpha, ...
    g);

savefig(sprintf('%s/lambda_from_ampl_and_d_%s.fig', home_figures_path, suffix));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/lambda_from_ampl_and_d_%s.pdf', home_figures_path, suffix));
