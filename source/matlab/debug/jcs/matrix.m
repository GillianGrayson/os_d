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
 
fn = sprintf('%s/exp_mtx_0_%s.txt', path, suffix);
data = importdata(fn);
exp_mtx_0 = zeros(sys_size);
for st_id_1 = 1:sys_size
    for st_id_2 = 1:sys_size
        index = (st_id_1 - 1) * sys_size + st_id_2;
        exp_mtx_0(st_id_1, st_id_2) = data(index, 1) + sqrt(-1) * data(index, 2);
    end
end

fn = sprintf('%s/exp_mtx_1_%s.txt', path, suffix);
data = importdata(fn);
exp_mtx_1 = zeros(sys_size);
for st_id_1 = 1:sys_size
    for st_id_2 = 1:sys_size
        index = (st_id_1 - 1) * sys_size + st_id_2;
        exp_mtx_1(st_id_1, st_id_2) = data(index, 1) + sqrt(-1) * data(index, 2);
    end
end

alfa = 5;
xi = 1;
gamma1 = 0.25 / alfa;
gamma2 = 0.25/alfa;
n_sr = 0;

a = zeros(N,N);
for n=1:(N-1)
    a(n,n+1) = sqrt(n);
end
a_ = a';

A1=sqrt(n_sr+1)*a;
A2=sqrt(n_sr)*a_;

ham = 1/2*(1/(alfa*alfa*alfa))*a_*a_*a*a;
ham_drv = (a_-a);
ham_diff = max(max(abs(ham - hamiltonian)))
ham_drv_diff = max(max(abs(ham_drv - hamiltonian_drv)))

H1=1/2*(1/(alfa*alfa*alfa))*a_*a_*a*a-(i/2)*(gamma1*A1'*A1+gamma2*A2'*A2);
H2=i*(a_-a);

f0 = 3.2;
F0 = f0*xi;
tau1 = 0.98*alfa;
T = 1.98*alfa;
t1 = tau1/xi;
t2 = T/xi - t1;

NNNt = 256;
h1 = t1/NNNt;
h2 = t2/NNNt;

NNNt1=NNNt;
NNNt2=NNNt;

ham_qj_0 = -i*(H1+F0*H2);
ham_qj_1 = -i*(H1);

G1 = expm(-i*(H1+F0*H2)*h1);
G2 = expm(-i*(H1)*h2);

G1_diff = max(max(abs(G1 - exp_mtx_0)))
G2_diff = max(max(abs(G2 - exp_mtx_1)))



    