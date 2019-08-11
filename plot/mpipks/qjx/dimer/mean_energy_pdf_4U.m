clear all;

home_figures_path = '/home/denysov/yusipov/os_d/figures';
data_path = '/data/biophys/denysov/yusipov/os_d/data/qjx';

sys_id = 0;
task_id = 1;
prop_id = 1;

ex_deep = 16;
rk_ns = 100000;
num_tp_periods = 101;
num_obs_periods = 100;

N = 500;

Us = [0.1; 0.1125; 0.125; 0.15];

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

num_trajectories = 1;
num_runs = 100;

sys_size = N + 1;

num_bins_x = 200;
num_bins_y = 200;

starts_x = [0.125; 0.5; 0.125; 0.5];
shifts_x = [0.37; 0.37; 0.37; 0.37]; 
starts_y = [0.53; 0.53; 0.125; 0.125];
shifts_y = [0.4; 0.4; 0.4; 0.4];

x_min = 0;
x_max = N;
x_shift = (x_max - x_min) / num_bins_x;
x_bins = linspace(x_min + 0.5 * x_shift, x_max - 0.5 * x_shift, num_bins_x) / N;


y_min = -2.0 * N;
y_max = 2.0 * N;
y_shift = (y_max - y_min) / num_bins_y;
y_bins = linspace(y_min + 0.5 * y_shift, y_max - 0.5 * y_shift, num_bins_y);

fig = figure;
for U_id = 1:size(Us, 1)
    
    dimer_prm_U = Us(U_id);
    
    mean_evo = zeros(num_obs_periods + 1, num_trajectories * num_runs);
    energy_evo = zeros(num_obs_periods + 1, num_trajectories * num_runs);
    
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
            dimer_prm_U, ...
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
            dimer_prm_U, ...
            dimer_prm_J, ...
            start_type, ...
            start_state);
        
        fn = sprintf('%s/mean_evo_%s.txt', path_to_folder, suffix);
        mean_evo_curr = importdata(fn);
        mean_evo(:, ss + 1: ss + num_trajectories) = mean_evo_curr;
        
        fn = sprintf('%s/energy_evo_%s.txt', path_to_folder, suffix);
        energy_evo_curr = importdata(fn);
        energy_evo(:, ss + 1: ss + num_trajectories) = energy_evo_curr;
    end
    
    z_data = zeros(num_bins_x, num_bins_y);
    
    total_num_trajectories = num_trajectories * num_runs;
    
    for tr_id = 1:total_num_trajectories
        
        xs = mean_evo(:, tr_id);
        ys = energy_evo(:, tr_id);
        
        for p_id = 1 : size(xs, 1)
            x = xs(p_id);
            y = ys(p_id);
            
            x_id = floor((x - x_min) / (x_max - x_min + eps) * num_bins_x) + 1;
            y_id = floor((y - y_min) / (y_max - y_min + eps) * num_bins_y) + 1;
            
            z_data(x_id, y_id) = z_data(x_id, y_id) + 1;
        end
        
    end
    
    z_data = z_data / (size(mean_evo, 1) * size(mean_evo, 2) * x_shift * y_shift);
    norm = sum(sum(z_data)) *  x_shift * y_shift
    z_data = z_data / max(max(z_data));
    
    subplot(2, 2, U_id);
    hLine = imagesc(x_bins, y_bins, z_data');
    
    set(gca, 'FontSize', 28);
    if U_id > 2
        xlabel('$n/N$', 'Interpreter', 'latex')
		set(gca,'xtick', [0.25 0.5 0.75 1])
    else
        xlabel('')
        set(gca,'xticklabel',{[]})
    end
    
    set(gca, 'FontSize', 28);
    if U_id == 1 || U_id == 3
        ylabel('$e$', 'Interpreter', 'latex');
		set(gca,'ytick', [-N 0 N])
    else
        ylabel('')
        set(gca,'yticklabel',{[]})
    end
    
    if U_id == 1
        colormap hot;
        h = colorbar;
        set(gca, 'FontSize', 28);
        title(h, '$PDF$', 'FontSize', 30, 'interpreter','latex');
        set(h, 'Position', [0.88 0.125 0.02 0.805]);
        set(h,'ytick', [0 0.25 0.5 0.75 1])
    else
        colorbar('off')
    end
    
    set(gca,'YDir','normal');
    set(gca, 'Position', [starts_x(U_id) starts_y(U_id) shifts_x(U_id) shifts_y(U_id)]);
    
    hold all;
   
end

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

savefig(sprintf('%s/mean_energy_pdf_4U_%s.fig', home_figures_path, suffix_save));

h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('%s/mean_energy_pdf_4U_%s.pdf', home_figures_path, suffix_save));

