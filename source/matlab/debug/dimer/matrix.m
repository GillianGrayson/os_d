clear all;

qj_seed = 0;
qj_mns = 1000000;

N = 10;
diss_type = 1;
diss_gamma = 0.1;
diss_phase = 0.0;
drv_type  = 0;
drv_ampl = 1.5;
drv_freq = 1.0;
drv_phase = 0.0;
start_type = 0;
start_state = 5;
prm_E = 1.0;
prm_U = 0.1;
prm_J = 1.0;

sys_size = N + 1;

path = "../../cpp/QJX/QJX";

suffix = sprintf("qjrnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%d_%0.4f_%0.4f_%0.4f)_start(%d_%d)_prm(%0.4f_%0.4f_%0.4f)", ...
    qj_seed, ...
    qj_mns, ...
    N, ...
    diss_type, ...
    diss_gamma, ...
    diss_phase, ...
    drv_type, ...
    drv_ampl, ...
    drv_freq, ...
    drv_phase, ...
    start_type, ...
    start_state, ...
    prm_E, ...
    prm_U, ...
    prm_J);

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
    