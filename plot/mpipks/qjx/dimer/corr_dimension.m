clear all;


home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'qjx_results';

data_path = sprintf('%s/%s', data_path, prefix);

sys_id = 0;
task_id = 2;
is_debug = 0;
is_pp = 1;
init_fn = '';
path = '';
num_threads = 1;
qj_deep = 4;
qj_num_tp_periods = 1000;
qj_num_obs_periods = 10;
qj_num_trajectories = 1;
qj_seed = 0;
qj_mns = 1000000;

type_lpn = 0;
eps_lpn = 1.0e-3;
delta_up_lpn = 1.0e-2;
delta_down_lpn = 1.0e-13;
is_obs_dump = 1;
is_adr_dump_sep = 0;
is_adr_dump_avg = 0;
is_evo_dump_sep = 1;
is_evo_dump_avg = 0;
dump_type = 0;
num_dumps = 1;
N = 200;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;
start_type = 0;
start_state = 0;
prm_E = 1.0;
prm_U = 0;
prm_J = 1.0;
cd_num_sub_steps = 25000; 
cd_dim = 1; 
cd_eps = 1e-08; 

sys_size = N + 1;

num_runs = 1;

U_begin = 0.1;
U_step = 0.6;
U_num = 2;

eps_start = 1.0e-8;
eps_per_decade = 10;
eps_num = 81;

m_begin = 1;
m_step = 1;
m_num = 5;

for U_id = 1:U_num
    
	U_id = U_id;
	U = U_begin + U_step * (U_id - 1)

	fig = figure;
	   
	for m_id = 1:m_num
        
        m = m_begin + m_step * (m_id - 1);
		
		epss = zeros(eps_num, 1);
		cis = zeros(eps_num, 1);
		
		for eps_id = 1:eps_num
			
			eps = eps_start * 10.0^((eps_id-1)/eps_per_decade);
			
			epss(eps_id) = eps;
            
            for run_id = 1:num_runs
                
                ss = (run_id - 1) * qj_num_trajectories;
                
                path_to_folder = sprintf('%s/main_%d_%d/qj_%d_%d_%d/ci_%d_%d_%0.10f/N_%d/diss_%d_%0.4f_%0.4f/drv_%d_%0.4f_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
                    data_path, ...
                    sys_id, ...
                    task_id, ...
                    qj_deep, ...
                    qj_num_tp_periods, ...
                    qj_num_obs_periods, ...
                    cd_num_sub_steps, ...
                    m, ...
                    eps, ...
                    N, ...
                    diss_type, ...
                    diss_gamma, ...
                    diss_phase, ...
                    drv_type, ...
                    drv_ampl, ...
                    drv_freq, ...
                    drv_phase, ...
                    prm_E, ...
                    U, ...
                    prm_J, ...
                    start_type, ...
                    start_state, ...
                    ss);
                
                suffix = sprintf('qjrnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
                    ss, ...
                    qj_mns, ...
                    N, ...
                    diss_type, ...
                    diss_gamma, ...
                    diss_phase, ...
                    drv_type, ...
                    drv_ampl, ...
                    drv_freq, ...
                    drv_phase, ...
                    prm_E, ...
                    U, ...
                    prm_J, ...
                    start_type, ...
                    start_state);
                
                
                fn = sprintf('%s/ci_%s.txt', path_to_folder, suffix);
                data = importdata(fn);
				
				cis(eps_id) = cis(eps_id) + data;
            end
            
            cis(eps_id) = cis(eps_id) / num_runs;
			
		end
		
		hLine = plot(log10(epss), log10(cis), 'LineWidth', 2);
		legend(hLine, sprintf('%d', m_id));
		
		set(gca, 'FontSize', 30);
		xlabel('$\log_{10}(\epsilon)$', 'Interpreter', 'latex');
		
		set(gca, 'FontSize', 30);
		ylabel('$\log_{10}(C)$', 'Interpreter', 'latex');
		
		hold all;
		
	end

	suffix = sprintf('qj(%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d).txt', ...
            num_runs * qj_num_trajectories, ...
            N, ...
            diss_type, ...
            diss_gamma, ...
            diss_phase, ...
            drv_type, ...
            drv_ampl, ...
            drv_freq, ...
            drv_phase, ...
            prm_E, ...
            U, ...
            prm_J, ...
            start_type, ...
            start_state);

	savefig(sprintf('%s/cd_%s.fig', home_figures_path, suffix));

end
