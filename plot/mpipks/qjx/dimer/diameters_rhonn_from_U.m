clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 0; 
task_id = 1;
prop_id = 1;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 100;
num_obs_periods = 100; 

N = 100;

diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;

dimer_drv_type = 1;
dimer_drv_ampl = 3.4;
dimer_drv_freq = 1.0;
dimer_drv_phase = 0.0;

dimer_prm_E = 0.0;
dimer_prm_J = 1.0;

start_type = 0;
start_state = 0;

ss = 0;
mns = 1000000;

num_trajectories = 10;
num_runs = 10;

sys_size = N + 1;

U_begin = 0.05;
U_step = 0.0025;
U_num = 101;

Us = zeros(U_num, 1);
states = linspace(1, sys_size, sys_size);

diameters = zeros(U_num, 1);

pks_lim = 1e-3;

for U_id = 1:U_num
	
	U_id = U_id;
	U = U_begin + U_step * (U_id - 1)
	Us(U_id) = U;
   
	curr_rho_avg = zeros(sys_size,1);
   
	for run_id = 1:num_runs
        
        ss = (run_id - 1) * num_trajectories;
        
        path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/N_%d/diss_%d_%0.4f_%0.4f/drv_%d_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
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
            dimer_drv_type, ...
            dimer_drv_ampl, ...
            dimer_drv_freq, ...
            dimer_drv_phase, ...
            dimer_prm_E, ...
            U, ...
            dimer_prm_J, ...
            start_type, ...
            start_state, ...
            ss);
        
        suffix = sprintf('rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
            ss, ...
            mns, ...
            N, ...
            diss_type, ...
            diss_gamma, ...
            diss_phase, ...
            dimer_drv_type, ...
            dimer_drv_ampl, ...
            dimer_drv_freq, ...
            dimer_drv_phase, ...
            dimer_prm_E, ...
            U, ...
            dimer_prm_J, ...
            start_type, ...
            start_state);
        
        path = sprintf('%s/adr_avg_%s.txt', path_to_folder, suffix);
        data = importdata(path);
        
        time_avg_rho = zeros(sys_size,1);
        for period_id = 1 : num_obs_periods + 1
            begin_index = (period_id - 1) * sys_size + 1;
            end_index = period_id * sys_size;
            curr_rho = data(begin_index : end_index);
            time_avg_rho = time_avg_rho + curr_rho;
        end
        time_avg_rho = time_avg_rho / (num_obs_periods + 1);
        
        curr_rho_avg = curr_rho_avg + time_avg_rho;
    end
    
    curr_rho_avg = curr_rho_avg / num_runs;
    
    norm = sum(curr_rho_avg);
    norm_diff = 1.0 - norm
	
	curr_rho_avg = smooth(curr_rho_avg, floor(N * 0.1));
    
    [pks, locs] = findpeaks(curr_rho_avg);
    num_peaks = size(pks, 1);
    
    del_ids = [];
    for i = 1:size(pks, 1)
        if pks(i) < pks_lim
            del_ids = vertcat(del_ids, i);
        end
    end
    pks(del_ids) = [];
    locs(del_ids) = [];
    
    if size(pks, 1) == 1 
        diameters(U_id) = 0;
    else
        diameters(U_id) = abs(locs(end) - locs(1)) / N;
    end
    
end



fig = figure;
hLine = plot(Us, diameters);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$D$', 'Interpreter', 'latex');
xlim([Us(1) Us(end)])

suffix_save = sprintf('N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_var_%0.4f)_start(%d_%d)', ...
        N, ...
        diss_type, ...
        diss_gamma, ...
        diss_phase, ...
        dimer_drv_type, ...
        dimer_drv_ampl, ...
        dimer_drv_freq, ...
        dimer_drv_phase, ...
        dimer_prm_E, ...
        dimer_prm_J, ...
        start_type, ...
        start_state);

savefig(sprintf('%s/diameters_rhonn_from_U_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/diameters_rhonn_from_U_%s.pdf', home_figures_path, suffix_save));

