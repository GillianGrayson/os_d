clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'qjx_results';

data_path = sprintf('%s/%s', data_path, prefix);

sys_id = 1;
task_id = 4;
prop_id = 0;
seed = 0;
mns = 1000000;
num_trajectories = 100;
num_tp_periods = 10000;
num_obs_periods = 10;
ex_deep = 10;
rk_ns = 10000;

N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0;
jcs_drv_part_1 = 0.98;
jcs_drv_part_2 = 1;
jcs_drv_ampl = 1.44;
jcs_prm_alpha = 5;
start_type = 0;
start_state = 0;

sys_size = N;

deep_num_steps = 128;

num_dumps = 2 * deep_num_steps * num_obs_periods + 1;
dump_periods = zeros(num_dumps, 1);

num_periods = (num_dumps - 1) / (2 * deep_num_steps);
dump_shift = 1 / (2 * deep_num_steps);
dump_periods(1) = 0;
for dump_id = 2:num_dumps
    dump_periods(dump_id) = dump_shift * (dump_id - 1);
end

num_runs = 100;

spec_begin = -3.0;
spec_end = 3.0;
spec_num = 500;
spec_shift = (spec_end - spec_begin) / spec_num;
specs = zeros(spec_num, 1);
for spec_id = 1 : spec_num
    specs(spec_id) = spec_begin + ((spec_id - 1) + 0.5) * spec_shift;
end

pdf_real = zeros(num_dumps, spec_num);
pdf_imag = zeros(num_dumps, spec_num);

non_inc_real = 0;
non_inc_imag = 0;

for run_id = 1:num_runs
    
    ss = (run_id - 1) * num_trajectories
    
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
        jcs_drv_part_1, ...
        jcs_drv_part_2, ...
        jcs_drv_ampl, ...
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
        jcs_drv_part_1, ...
        jcs_drv_part_2, ...
        jcs_drv_ampl, ...
        jcs_prm_alpha, ...
        start_type, ...
        start_state);
    
    
    path = sprintf('%s/spec_evo_%s.txt', path_to_folder, suffix);
    data = importdata(path);
    data_size = size(data, 1);
    
    for dump_id = 1:num_dumps
	
		pdf_real_curr = zeros(spec_num, 1);
		pdf_imag_curr = zeros(spec_num, 1);
        
        for tr_id = 1:num_trajectories
            
            col_id_re = (tr_id - 1) * 2 + 1;
            col_id_im = (tr_id - 1) * 2 + 2;
            
            curr_real = data(dump_id, col_id_re);
            if ((curr_real >= spec_begin) && (curr_real <= spec_end))
                bin_id = floor((curr_real - spec_begin) / spec_shift + 1.0e-10) + 1;
                pdf_real_curr(bin_id) = pdf_real_curr(bin_id) + 1;
            else
                non_inc_real = non_inc_real + 1;
            end
            
            curr_imag = data(dump_id, col_id_im);
            if ((curr_imag >= spec_begin) && (curr_imag <= spec_end))
                bin_id = floor((curr_imag - spec_begin) / spec_shift + 1.0e-10) + 1;
                pdf_imag_curr(bin_id) = pdf_imag_curr(bin_id) + 1;
            else
                non_inc_real = non_inc_real + 1;
            end
            
        end
        
        real_count = sum(pdf_real_curr);
        pdf_real_curr = pdf_real_curr / (real_count * spec_shift);
        real_norm = sum(pdf_real_curr) * spec_shift;
        pdf_real_curr = pdf_real_curr / max(pdf_real_curr);
        pdf_real(dump_id, :) = pdf_real_curr;
        
        imag_count = sum(pdf_imag_curr);
        pdf_imag_curr = pdf_imag_curr / (imag_count * spec_shift);
        imag_norm = sum(pdf_imag_curr) * spec_shift;
        pdf_imag_curr = pdf_imag_curr / max(pdf_imag_curr);
        pdf_imag(dump_id, :) = pdf_imag_curr;
        
    end
    
end

non_inc_real = non_inc_real
non_inc_imag = non_inc_imag

pdf_real = pdf_real';
pdf_imag = pdf_imag';

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
    jcs_drv_part_1, ...
    jcs_drv_part_2, ...
	jcs_drv_ampl, ...
    jcs_prm_alpha, ...
    start_type, ...
    start_state);

fig = figure;
hLine = imagesc(dump_periods, specs, pdf_real);
set(gca, 'FontSize', 30);
xlabel('$f_0$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$Re(\xi)$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF');
set(gca,'YDir','normal');

savefig(sprintf('%s/deep_real_spec_from_ampl_%s.fig', home_figures_path, suffix));

close(fig);

fig = figure;
hLine = imagesc(dump_periods, specs, pdf_imag);
set(gca, 'FontSize', 30);
xlabel('$f_0$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$Im(\xi)$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'PDF');
set(gca,'YDir','normal');

savefig(sprintf('%s/deep_imag_spec_from_ampl_%s.fig', home_figures_path, suffix));

close(fig);