clear all;

tr_id = 1;
num_sub_steps = 1024;

sys_id = 3;
task_id = 4;
prop_id = 0;

seed = 1;
mns = 1000000;

Nc = 8;

W_seed = 10;
W_mns = 1000000;

diss_type = 1;
diss_phase = 0.0;
diss_gamma = 0.1;

prm_W = 2.0;
prm_U = 10.0;
prm_J = 1.0;

start_type = 0;
start_state = 49;

delta_lpn_s = 1e-5;
delta_lpn_f_h = 1000000;
delta_lpn_f_l = 1e-12;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('setup(%d_%d_%d)_rnd(%d_%d)_Nc(%d)_rnd(%d_%d)_diss(%d_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
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
    start_state);

fn = sprintf('%s/periods_%s.txt', data_path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

num_periods = (num_dumps - 1) / (1 * num_sub_steps);
dump_shift = 1 / (1 * num_sub_steps);
dump_periods(1) = 0;
for dump_id = 2:num_dumps
    dump_periods(dump_id) = dump_shift * (dump_id - 1);
end

N = size(dump_periods, 1);
fft_bins = ((1/N:1/N:1) * num_sub_steps).';

fn = sprintf('%s/spec_evo_%s.txt', data_path, suffix);
spec_evo_data = importdata(fn);

fn = sprintf('%s/imbalance_evo_%s.txt', data_path, suffix);
imbalance_evo_data = importdata(fn);

spec_evo = spec_evo_data(:, 2 * (tr_id - 1) + 1);
imbalance_evo = imbalance_evo_data(:, tr_id);

spec_fft = fft(spec_evo);
imbalance_fft = fft(imbalance_evo);


fig = figure;
hLine = plot(fft_bins, abs(spec_fft));
xlabel('frequency (Hz)', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('magnitude', 'Interpreter', 'latex');
xlim([fft_bins(1) fft_bins(end)])
set(gca, 'YScale', 'log')
hold all;
hLine = plot(fft_bins, abs(imbalance_fft));

legends = cell(2,1);
legends{1} = sprintf('energy W=%0.2f U=%0.2f', prm_W, prm_U);
legends{2} = sprintf('imbalance W=%0.2f U=%0.2f', prm_W, prm_U);
legend(legends, 'FontSize', 14, 'Location', 'North');

suffix = sprintf('ns(%d)_prm(%0.4f_%0.4f_%0.4f)', ...
    Nc, ...
    prm_W, ...
    prm_U, ...
    prm_J);

savefig(sprintf('tr_fft_%s.fig', suffix));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', sprintf('tr_fft_%s.pdf', suffix));

propertyeditor(fig)



