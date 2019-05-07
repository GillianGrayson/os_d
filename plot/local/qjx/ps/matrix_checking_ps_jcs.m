clear all;

seed = 1;
mns = 1000000;
num_trajectories = 1;

N = 10;
diss_type = 0;
diss_gamma = 0.1;
diss_phase = 0.0;
jcs_drv_part_1 = 0.98;
jcs_drv_part_2 = 1.0;
jcs_drv_ampl = 3.2;
jcs_prm_alpha = 5.0;

ps_num_spins = 1;
ps_num_spins_states = 2^ps_num_spins;
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

%matrix_name = "hamiltonian";
%matrix_name = "hamiltonian_drv";
%matrix_name = "dissipator_0";
%matrix_name = "hamiltonians_qj_0";
matrix_name = "hamiltonians_qj_1";

num_bins = 201;

path = "../../../../source/cpp/QJX/QJX";

suffix = sprintf("rnd(%d_%d)_N(%d)_diss(%d_%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f)_start(%d_%d)", ...
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

fn = sprintf('%s/%s_%s.txt', path, matrix_name, suffix);
data = importdata(fn);

jcs_mtx = zeros(N);
for i = 1:N
    for j = 1:N
        index = (i-1) * N + j;
        if size(data, 2) == 1
            jcs_mtx(i, j) = data(index);
        else
            jcs_mtx(i, j) = complex(data(index, 1), data(index, 2));
        end
    end
end

sigma_0 = eye(ps_num_spins_states);

jcs_mtx_kron = kron(sigma_0, jcs_mtx);

suffix = sprintf("rnd(%d_%d)_s(%d)_nps(%d)_diss(%0.4f_%0.4f)_drv(%0.4f_%0.4f_%0.4f)_prm(%0.4f_%0.4f_%0.4f)_start(%d_%d)", ...
    seed, ...
    mns, ...
    ps_num_spins, ...
    ps_num_photons_states, ...
    ps_diss_w, ...
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

diff = max(max(abs(jcs_mtx_kron - ps_mtx)))




