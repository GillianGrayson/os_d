clear all;
addpath('../../../../../os_lnd/source/matlab/lib')

figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 3;
task_id = 7;
prop_id = 0;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 500;
num_obs_periods = 500;

T = 10;

Nc = 8;

num_rnd_obs = 1;
rnd_obs_seed = 100;
rnd_obs_type = 2;

W_seed_s = 4;
W_seed_num = 1;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

U = 1;
J = 1;

start_type = 0;
start_state = 0;

lpn_type = -1;
lpn_delta_f_h = 1e-6;
lpn_delta_s = 1e-6;
lpn_delta_f_l = 1e-6;

ss = 0;
mns = 1000000;

num_trajectories = 10000;
num_target_trajectories = num_trajectories / 2;
num_runs = 1;

W_begin = 20.0;
W_step = 0.2;
W_num = 1;
Ws = zeros(W_num, 1);

all_lambdas = zeros(W_num, W_seed_num * num_runs * num_target_trajectories);

for W_id = 1:W_num
    
    W = W_begin + W_step * (W_id - 1)
    Ws(W_id) = W;
       
    for W_seed_id = 1 : W_seed_num
        
        W_seed = W_seed_s + W_seed_id - 1;
        
        for run_id = 1 : num_runs
            
            ss = (run_id - 1) * num_trajectories;
            
            path_to_folder = sprintf('%s/main_%d_%d_%d/T_%0.4f/lpn_%d_%0.4f_%0.4f_%0.4f/run_%d_%d_%d_%d/Nс_%d/obs_%d_%d_%d/rnd_%d_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/ss_%d', ...
				data_path, ...
				sys_id, ...
				task_id, ...
				prop_id, ...
				T, ...
				lpn_type, ...
				log10(lpn_delta_s), ...
				log10(lpn_delta_f_h), ...
				log10(lpn_delta_f_l), ...
				ex_deep, ...
				rk_ns, ...
				num_tp_periods, ...
				num_obs_periods, ...
				Nc, ...
				num_rnd_obs, ...
				rnd_obs_seed, ...
				rnd_obs_type, ...
				W_seed, ...
				W_mns, ...
				diss_type, ...
				diss_phase, ...
				diss_gamma, ...
				W, ...
				U, ...
				J, ...
				start_type, ...
				start_state, ...
				ss);
            
            suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
                sys_id, ...
                task_id, ...
                prop_id, ...
                ss, ...
                mns, ...
                Nc, ...
                W_seed, ...
                W_mns, ...
                diss_type, ...
                diss_phase, ...
                diss_gamma, ...
                W, ...
                U, ...
                J, ...
                start_type, ...
                start_state, ...
                lpn_type, ...
                log10(lpn_delta_s), ...
                log10(lpn_delta_f_h), ...
                log10(lpn_delta_f_l));
            s_index = (W_seed_id - 1) * num_runs * num_target_trajectories + (num_runs - 1) * num_target_trajectories + 1;
            f_index = (W_seed_id - 1) * num_runs * num_target_trajectories + (num_runs - 1) * num_target_trajectories + num_target_trajectories;
            
            fn = sprintf('%s/lambda_%s.txt', path_to_folder, suffix);
            curr_data = importdata(fn);
            curr_lambdas = curr_data(num_target_trajectories + 1 : end);
            all_lambdas(W_id, s_index : f_index) = curr_lambdas;
            
            fn = sprintf('%s/num_renorms_%s.txt', path_to_folder, suffix);
            curr_data = importdata(fn);
            curr_renorms = curr_data(num_target_trajectories + 1 : end);
            all_renorms(W_id, s_index : f_index) = curr_renorms;
        end 
    end
    
    pdf.x_num_bins = 200;
    pdf.x_label = '$\lambda$';
    
    lambdas = all_lambdas(W_id, :)';
	ololo = size(lambdas)
    
    suffix_save =  sprintf('T(%0.4f)_setup(%d_%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
        T, ...
        sys_id, ...
        task_id, ...
        prop_id, ...
        Nc, ...
        W_seed, ...
        W_mns, ...
        diss_type, ...
        diss_phase, ...
        diss_gamma, ...
        W, ...
        U, ...
        J, ...
        start_type, ...
        start_state, ...
        lpn_type, ...
        log10(lpn_delta_s), ...
        log10(lpn_delta_f_h), ...
        log10(lpn_delta_f_l));
    
    pdf.x_bin_s = min(lambdas) - 1e-6;
    pdf.x_bin_f = max(lambdas) + 1e-6;
    pdf = oqs_pdf_1d_setup(pdf);
    pdf = oqs_pdf_1d_update(pdf, lambdas);
    pdf = oqs_pdf_1d_release(pdf);
    fig = oqs_pdf_1d_plot(pdf);
    fn_fig = sprintf('%s/lambdas_%s', figures_path, suffix_save);
    oqs_save_fig(fig, fn_fig);
    
    fn_txt = sprintf('%s/lambdas_%s.txt', figures_path, suffix_save);
    fid = fopen(fn_txt,'wt');
    for x_id = 1:size(lambdas, 1)
        fprintf(fid,'%0.16e\n', lambdas(x_id));
    end
    fclose(fid);
    
end
