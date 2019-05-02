clear all;

sys_id = 0;
task_id = 5;
prop_id = 0;

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

num_tr = 1;

is_mean = 1;

%data_path = '../../../../source/cpp/QJX/QJX';
data_path = 'C:\Users\user\Desktop\New folder\U(0.50)\T4_ds(1e-6)'

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

states = linspace(1, size_sys, size_sys) / size_sys;

fn = sprintf('%s/lambda_%s.txt', data_path, suffix);
lambda_data = importdata(fn);

mean_l = mean(lambda_data(2:end))
