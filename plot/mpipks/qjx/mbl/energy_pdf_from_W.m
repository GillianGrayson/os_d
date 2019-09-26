clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 3; 
task_id = 7;
prop_id = 0;

ex_deep = 16;
rk_ns = 10000;
num_tp_periods = 1000;
num_obs_periods = 1000;

Nc = 8;

W_seed = 1;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

U = 1;
J = 1;

start_type = 0;
start_state = 49;

lpn_delta_f = 1000;
lpn_delta_s = 0.00001;

ss = 0;
mns = 1000000;

num_trajectories = 1000;
num_target_trajectories = num_trajectories / 2;
num_runs = 1;

spec_min = -10;
spec_max = 10;
spec_num_bins = 400;
spec_shift = (spec_max - spec_min) / spec_num_bins;
spec_bins = linspace(spec_min + 0.5 * spec_shift, spec_max - 0.5 * spec_shift, spec_num_bins);

W_begin = 0.2;
W_step = 0.2;
W_num = 100;
Ws = zeros(W_num, 1);

spec_2d_pdf = zeros(W_num, spec_num_bins);

for W_id = 1:W_num
	
	W = W_begin + W_step * (W_id - 1)
	Ws(W_id) = W;
    
    spec_pdf = zeros(spec_num_bins, 1);
    non_inc = 0;

    for run_id = 1:num_runs

        ss = (run_id - 1) * num_trajectories;

        path_to_folder = sprintf('%s/main_%d_%d_%d/run_%d_%d_%d_%d/Nс_%d/rnd_%d_%d/diss_%d_%0.4f_%0.4f/prm_%0.4f_%0.4f_%0.4f/start_%d_%d/lpn_%0.4f_%0.4f/ss_%d', ...
            data_path, ...
            sys_id, ...
            task_id, ...
            prop_id, ...
            ex_deep, ...
            rk_ns, ...
            num_tp_periods, ...
            num_obs_periods, ...
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
            log10(lpn_delta_s), ...
            log10(lpn_delta_f), ...
            ss);

        suffix = sprintf('rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%0.4f_%0.4f)', ...
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
            log10(lpn_delta_s), ...
            log10(lpn_delta_f));

        fn = sprintf('%s/spec_evo_%s.txt', path_to_folder, suffix);
        spec_data = importdata(fn);
        
        for tr_id = 1:num_target_trajectories
            for p_id = 1:size(spec_data, 1)
                
				curr_spec = spec_data(p_id, 2*(tr_id-1) + 1);
				
                if curr_spec < spec_max && curr_spec >= spec_min
                    spec_id = floor((curr_spec - spec_min) / (spec_max - spec_min + eps) * spec_num_bins) + 1;
                    spec_pdf(spec_id) = spec_pdf(spec_id) + 1;
                else
                    non_inc = non_inc + 1;
                end
            end
        end
    end
    
    non_inc = non_inc
    
    sum_spec_pdf = sum(spec_pdf);
    spec_pdf = spec_pdf / (sum_spec_pdf * spec_shift);
    norm = sum(spec_pdf) * spec_shift;
    norm_diff = 1.0 - norm
    
    spec_2d_pdf(W_id, :) = spec_pdf;
end

for W_id = 1:W_num
    curr_max = max(spec_2d_pdf(W_id, :));
    spec_2d_pdf(W_id, :) = spec_2d_pdf(W_id, :) / curr_max;
end

fig = figure;
hLine = imagesc(Ws, spec_bins, spec_2d_pdf');
set(gca, 'FontSize', 30);
xlabel('$W$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('energy', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '');
set(gca,'YDir','normal');

suffix_save =  sprintf('Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(var_%0.4f_%0.4f)_start(%d_%d)_lpn(%0.4f_%0.4f)', ...
    Nc, ...
    W_seed, ...
    W_mns, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    U, ...
    J, ...
    start_type, ...
    start_state, ...
    log10(lpn_delta_s), ...
    log10(lpn_delta_f));

savefig(sprintf('%s/energy_pdf_from_W_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/energy_pdf_from_W_%s.pdf', home_figures_path, suffix_save));

