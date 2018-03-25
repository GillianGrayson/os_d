clear all;

sys_id = 1;
task_id = 1;
prop_id = 0;

seed = 0;
mns = 1000000;

N = 10;

diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;

jcs_drv_part_1 = 0.98;
jcs_drv_part_2 = 1.0;
jcs_drv_ampl = 3.2;
jcs_prm_alpha = 5.0;

start_type = 0;
start_state = 0;

sys_size = N;

path = "../../../cpp/QJX/QJX";

suffix = sprintf("config(%d_%d_%d)_rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)", ...
    sys_id, ...
    task_id, ...
    prop_id, ...
    seed, ...
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

fn = sprintf('%s/dissipator_drv_%s.txt', path, suffix);
data = importdata(fn);
dissipator = zeros(sys_size);
for st_id_1 = 1:sys_size
    for st_id_2 = 1:sys_size
        index = (st_id_1 - 1) * sys_size + st_id_2;
        dissipator(st_id_1, st_id_2) = data(index, 1) + sqrt(-1) * data(index, 2);
    end
end

fn = sprintf('%s/hamiltonian_%s.txt', path, suffix);
data = importdata(fn);
hamiltonian = zeros(sys_size);
for st_id_1 = 1:sys_size
    for st_id_2 = 1:sys_size
        index = (st_id_1 - 1) * sys_size + st_id_2;
        hamiltonian(st_id_1, st_id_2) = data(index, 1);
    end
end

fn = sprintf('%s/hamiltonian_drv_%s.txt', path, suffix);
data = importdata(fn);
hamiltonian_drv = zeros(sys_size);
for st_id_1 = 1:sys_size
    for st_id_2 = 1:sys_size
        index = (st_id_1 - 1) * sys_size + st_id_2;
        hamiltonian_drv(st_id_1, st_id_2) = data(index, 1);
    end
end

fn = sprintf('%s/hamiltonians_qj_0_%s.txt', path, suffix);
data = importdata(fn);
hamiltonians_qj_0 = zeros(sys_size);
for st_id_1 = 1:sys_size
    for st_id_2 = 1:sys_size
        index = (st_id_1 - 1) * sys_size + st_id_2;
        hamiltonians_qj_0(st_id_1, st_id_2) = data(index, 1) + sqrt(-1) * data(index, 2);
    end
end
    
fn = sprintf('%s/hamiltonians_qj_1_%s.txt', path, suffix);
data = importdata(fn);
hamiltonians_qj_1 = zeros(sys_size);
for st_id_1 = 1:sys_size
    for st_id_2 = 1:sys_size
        index = (st_id_1 - 1) * sys_size + st_id_2;
        hamiltonians_qj_1(st_id_1, st_id_2) = data(index, 1) + sqrt(-1) * data(index, 2);
    end
end
    