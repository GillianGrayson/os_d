clear all;

sys_id = 0;
task_id = 4;
prop_id = 0;

seed = 0;
mns = 1000000;
N = 100;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;
prm_E = 1.0;
prm_U = 0.01;
prm_J = 1.0;
start_type = 0;
start_state = 0;

T = 2 * pi / drv_freq;

cd_dump_deep = 1;
cd_num_sub_steps = 1000;

size_sys = N + 1;

tr_id = 0;

epsilon = 1e-8;
dim = 1;

data_path = '../../../../source/cpp/QJX/QJX';

suffix = sprintf('config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)', ...
    sys_id, ...
    task_id, ...
    prop_id, ...
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
dump_periods = dump_periods(2:end);
num_dumps = size(dump_periods, 1);

num_periods = (num_dumps - 1) / (2 * cd_num_sub_steps);
dump_shift = 1 / (2 * cd_num_sub_steps);
dump_periods(1) = 0;
for dump_id = 2:num_dumps
    dump_periods(dump_id) = dump_shift * (dump_id - 1);
end
    
adr = zeros(size_sys, num_dumps);

states = linspace(1, size_sys, size_sys) / size_sys;

fn = sprintf('%s/mean_evo_%s.txt', data_path, suffix);
mean_evo_data = importdata(fn);
mean_evo_data = mean_evo_data(2:end);
mean_evo = mean_evo_data(:, tr_id + 1);


num_points = num_dumps - dim + 1;
integral = 0;
for point_id_1 = 1:num_points

    for point_id_2 = 1:num_points
        
        vec_1 = mean_evo(point_id_1 : point_id_1 + dim - 1);
        vec_2 = mean_evo(point_id_2 : point_id_2 + dim - 1);
        
        if (point_id_1 ~= point_id_2)
            for dim_id = 1:dim
                curr_diff_vector = vec_1 - vec_2;            
            end
            
            curr_norm = norm(curr_diff_vector) / size_sys;
            
            if curr_norm < epsilon
                
                point_id_1 = point_id_1
                point_id_2 = point_id_2
                
                integral = integral + 1;
            end
        end
    end
end

integral = integral / (num_points * (num_points - 1));

