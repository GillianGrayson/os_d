clear all;

seed = 0;
mns = 1000000;

ps_num_spins = 1;
ps_num_photons_states = 10;
ps_drv_part_1 = 0.98; 
ps_drv_part_2 = 1.00; 
ps_drv_ampl = 3.2;
ps_prm_alpha = 5;
ps_prm_d = 0.0;
ps_prm_g = 0.0;
ps_diss_w = 0.05;
start_type = 0;
start_state = 0;

ps_num_spins_states = 2^ps_num_spins;
sys_size = ps_num_photons_states * ps_num_spins_states;

sigma_0 = [1 0; 0 1];
sigma_x = [0 1; 1 0];
sigma_y = [0 -i; i 0];
sigma_z = [1 0; 0 -1];
sigma_minus = [0 0; 1 0];
sigma_plus = [0 1; 0 0];

s_eye = eye(ps_num_spins_states);
p_eye = eye(ps_num_photons_states);

a_std_p = zeros(ps_num_photons_states);
for n = 1 : (ps_num_photons_states-1)
    a_std_p(n, n+1) = sqrt(n);
end
a_dag_p = a_std_p';

diss{1} = kron(s_eye, a_std_p);

J_minus = sigma_minus * 0.5;
J_plus = sigma_plus * 0.5;

H_1_p = 1 / (2 * ps_prm_alpha^3) * a_dag_p * a_dag_p * a_std_p * a_std_p;
H_1 = kron(s_eye, H_1_p);

H2 = kron(s_eye, a_dag_p) * kron(J_minus, p_eye) + kron(J_plus, p_eye) * kron(s_eye, a_std_p);



A1=sqrt(n_sr+1)*a;
A2=sqrt(n_sr)*a_;
%??????????? ???????????? H(t)=H1+F(t)*H2
H1=-(i/2)*(gamma1*A1'*A1+gamma2*A2'*A2);
H2=i*(a_-a);

%matrix_name = "hamiltonian";
%matrix_name = "hamiltonian_drv";
matrix_name = "dissipator_1";
%matrix_name = "hamiltonians_qj_0";
%matrix_name = "hamiltonians_qj_1";


path = "../../../../source/cpp/QJX/QJX";

suffix = sprintf("rnd(%d_%d)_s(%d)_nps(%d)_diss(1_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
    seed, ...
    mns, ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    ps_diss_w, ...
    ps_drv_part_1, ...
    ps_drv_part_2, ...
    ps_drv_ampl, ...
    ps_prm_alpha, ...
    ps_prm_d, ...
    ps_prm_g, ...
    start_type, ...
    start_state);

fn = sprintf('%s/%s_%s.txt', path, matrix_name, suffix);
data = importdata(fn);

ps_N = ps_num_photons_states * ps_num_spins_states;
ps_mtx = zeros(ps_N);
for i = 1:ps_N
    for j = 1:ps_N
        index = (i-1) * ps_N + j;
        if size(data, 2) == 1
            ps_mtx(i, j) = data(index);
        else
            ps_mtx(i, j) = complex(data(index, 1), data(index, 2));
        end
    end
end

sigma_minus = [0 0; 1 0];
sigma_0 = eye(ps_num_photons_states);

res = kron(sigma_minus, sigma_0);

diff = max(max(abs(res - ps_mtx)))




