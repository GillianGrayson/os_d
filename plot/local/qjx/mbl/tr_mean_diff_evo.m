clear all;

seed = 1;
mns = 1000000;

Nc = 8;

W_seed = 10;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

prm_W = 0.2;
prm_U = 1.0;
prm_J = 1.0;

start_type = 0;
start_state = 49;

delta_lpn_s = 1e-6;
delta_lpn_f = 1e-3;

deep_dump = 0;
deep_num_steps = 100;

num_plot_traj = 50;

fig = figure;

num_trajectories = 50;


data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)_lpn(%0.4f_%0.4f)', ...
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
    log10(delta_lpn_s), ...
    log10(delta_lpn_f));

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

fn = sprintf('%s/spec_evo_%s.txt', data_path, suffix);
spec_evo_data = importdata(fn);

for tr_id = 1 : num_plot_traj
    
    spec_evo_base = spec_evo_data(:, 2 * (tr_id - 1) + 1);
    spec_evo_var = spec_evo_data(:, 2 * (tr_id - 1 + num_trajectories) + 1);

    spec_evo_diff = abs(spec_evo_base - spec_evo_var);

    hLine = plot(dump_periods, spec_evo_diff);
    legend(hLine, sprintf('tr=%d', tr_id))
    set(gca, 'FontSize', 30);
    xlabel('$t/T$', 'Interpreter', 'latex');
    xlim([dump_periods(1) dump_periods(end)])
    set(gca, 'FontSize', 30);
    ylabel('$|\Delta|$', 'Interpreter', 'latex');
    hold all;
    set(gca, 'YScale', 'log')
    ylim([1e-12, 10]);
    
end

set(gca, 'YScale', 'log')
ylim([1e-12, 10]);
plot([dump_periods(1) dump_periods(end)], [delta_lpn_f delta_lpn_f], 'LineStyle', '--', 'color', 'black');
plot([dump_periods(1) dump_periods(end)], [delta_lpn_s delta_lpn_s], 'LineStyle', '--', 'color', 'black');

fn = sprintf('%s/lambda_%s.txt', data_path, suffix);
lambda_data = importdata(fn);
mean_lambda = mean(lambda_data(num_trajectories + 1:end))

propertyeditor(fig)

