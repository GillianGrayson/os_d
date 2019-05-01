clear all;

seed = 1;
mns = 1000000;
N = 200;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;
prm_E = 1.0;
prm_U = 0.5;
prm_J = 1.0;
start_type = 0;
start_state = 0;

T = 2 * pi / drv_freq;

cd_dump_deep = 1;
cd_num_sub_steps = 100;

size_sys = N + 1;

num_tr = 7;

is_mean = 1;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
    seed, ...
    mns, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase, ...
    prm_E, ...
    prm_U, ...
    prm_J, ...
    start_type, ...
    start_state);

fn = sprintf('%s/periods_%s.txt', data_path, suffix);
dump_periods = importdata(fn);
num_dumps = size(dump_periods, 1);

if(cd_dump_deep == 1)
    num_periods = (num_dumps - 1) / (2 * cd_num_sub_steps);
    dump_shift = 1 / (2 * cd_num_sub_steps);
    dump_periods(1) = 0;
    for dump_id = 2:num_dumps
        dump_periods(dump_id) = dump_shift * (dump_id - 1);
    end
end


fn = sprintf('%s/lambda_evo_%s.txt', data_path, suffix);
evo_data = importdata(fn);



fig = figure;
for tr_id = 1:num_tr

    evo_curr = evo_data(:, tr_id + 1);
    hLine = plot(dump_periods, evo_curr);
    set(gca, 'FontSize', 30);
    xlabel('$t/T$', 'Interpreter', 'latex');
    xlim([1 dump_periods(end)])
    ylim([-1 1])
    set(gca, 'XScale', 'log')
    set(gca, 'FontSize', 30);
    ylabel('$\lambda$', 'Interpreter', 'latex');
    hold all;
    
end

propertyeditor(fig)

fn = sprintf('%s/num_renorms_%s.txt', data_path, suffix);
num_renorms_data = importdata(fn);
num_renorms = sum(num_renorms_data(2:num_tr+1)) / num_tr


