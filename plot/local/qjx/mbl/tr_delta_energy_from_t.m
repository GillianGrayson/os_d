clear all;

T = 0.1;

sys_id = 3;
task_id = 8;
prop_id = 0;

is_jump = 1;

seed = 1;
mns = 1000000;

Nc = 8;

W_seed = 10;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

prm_W = 20.0;
prm_U = 1.0;
prm_J = 1.0;

start_type = 0;
start_state = 0;

lpn_type = -1;
delta_lpn_s = 1e-6;
delta_lpn_f_h = 10;
delta_lpn_f_l = 1e-16;

deep_dump = 1;
deep_num_steps = 100;

num_target_traj = 10;
num_plot_traj = 10;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%d_%0.4f_%0.4f_%0.4f)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    seed, ...
    mns, ...
    Nc, ...
    W_seed, ...
    W_mns, ...
    diss_type, ...
    diss_phase, ...
    diss_gamma, ...
    prm_W, ...
    prm_U, ...
    prm_J, ...
    start_type, ...
    start_state, ...
    lpn_type, ...
    log10(delta_lpn_s), ...
    log10(delta_lpn_f_h), ...
    log10(delta_lpn_f_l));

fn = sprintf('%s/periods_%s.txt', data_path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

if(deep_dump == 1)
    num_periods = (num_dumps - 1) / (1 * deep_num_steps);
    dump_shift = 1 / (1 * deep_num_steps);
    dump_periods(1) = 0;
    for dump_id = 2:num_dumps
        dump_periods(dump_id) = dump_shift * (dump_id - 1);
    end
end

dump_periods = dump_periods * T;

fn = sprintf('%s/spec_lpn_evo_%s.txt', data_path, suffix);
spec_evo_data = importdata(fn);

fn = sprintf('%s/lambda_%s.txt', data_path, suffix);
lambda_data = importdata(fn);

fig = figure;
for tr_id = 1 : num_plot_traj
    
    spec_evo_base = spec_evo_data(:, 2 * (tr_id - 1) + 1);
    spec_evo_var = spec_evo_data(:, 2 * (tr_id - 1 + num_target_traj) + 1);

    spec_evo_diff = abs(spec_evo_base - spec_evo_var);

    hLine = plot(dump_periods, spec_evo_diff);
    legend(hLine, sprintf('\\lambda=%0.2e', lambda_data(num_target_traj + tr_id)), 'FontSize', 14)
    set(gca, 'FontSize', 30);
    xlabel('$t/T$', 'Interpreter', 'latex');
    xlim([dump_periods(1) dump_periods(end)])
    set(gca, 'FontSize', 30);
    ylabel('$|\Delta(\varepsilon)|$', 'Interpreter', 'latex');
    hold all;
    set(gca, 'YScale', 'log')
    ylim([1e-12, 10]);
    
    
    if (is_jump == 1)
        color = get(hLine,'Color');
        
        fn = sprintf('%s/jump_times_%d_%s.txt', data_path, (tr_id - 1), suffix);
        jump_data_1 = importdata(fn);
        for j_id = 1:size(jump_data_1, 1)
            plot([jump_data_1(j_id) jump_data_1(j_id)], [1e-20 1e+1], 'Color', color)
        end
        
        fn = sprintf('%s/jump_times_%d_%s.txt', data_path, (tr_id - 1 + num_target_traj), suffix);
        jump_data_2 = importdata(fn);
        for j_id = 1:size(jump_data_2, 1)
            plot([jump_data_2(j_id) jump_data_2(j_id)], [1e-20 1e+1], 'Color', color)
        end
    end
    
end

mean_lambda = mean(lambda_data(num_target_traj + 1:end));

set(gca, 'YScale', 'log')
h = plot([dump_periods(1) dump_periods(end)], [delta_lpn_f_h delta_lpn_f_h], 'LineStyle', '--', 'color', 'black');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h = plot([dump_periods(1) dump_periods(end)], [delta_lpn_f_l delta_lpn_f_l], 'LineStyle', '--', 'color', 'black');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h = plot([dump_periods(1) dump_periods(end)], [delta_lpn_s delta_lpn_s], 'LineStyle', '--', 'color', 'black');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([delta_lpn_f_l * 0.1, delta_lpn_f_h * 10]);

title_str = sprintf('$ns=%d$ \\quad $W=%d$ \\quad $U=%d$ \\quad $J=%d$ \\quad $log_{10}\\Delta_s=%.2f$ \\quad $log_{10}\\Delta^h_f=%.2f$ \\quad $log_{10}\\Delta^l_f=%.2f$ \\quad $\\lambda=%0.2e$', ...
    Nc, ...
    prm_W, ...
    prm_U, ...
    prm_J, ...
    log10(delta_lpn_s), ...
    log10(delta_lpn_f_h), ...
    log10(delta_lpn_f_l), ...
    mean_lambda);
title(title_str, 'interpreter', 'latex', 'FontSize', 20)

propertyeditor(fig)
b = gca; legend(b,'off');

suffix = sprintf('ns(%d)_prm(%0.4f_%0.4f_%0.4f)_lpn(%d_%0.4f_%0.4f)', ...
    Nc, ...
    prm_W, ...
    prm_U, ...
    prm_J, ...
    0, ...
    log10(delta_lpn_s), ...
    log10(delta_lpn_f_h));

savefig(sprintf('tr_delta_energy_evo_%s.fig', suffix));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('tr_delta_energy_evo_%s.pdf', suffix));



