clear all;

delta = 0.1;
omega = 1.0;

amplitude = 0.6;
freq = 0.07;

gamma = 0.01;

% photonic mode subsystem
num_photons = 1;
photonic_size = num_photons + 1;
a_std = zeros(photonic_size, photonic_size);
for n = 1 : photonic_size - 1
    a_std(n, n+1) = sqrt(n);
end
a_dag = a_std';

% system size (we have only one spin)
N = 2 * photonic_size;
N2 = N * N;

% Matrixes
s_0 = eye(2);
s_x = [0 1; 1 0];
s_y = [0 -1i; 1i 0];
s_z = [1 0; 0 -1];
s_minus = [0 0; 1 0];
s_plus = [0 1; 0 0];

% Hamiltonian
sigma_z = kron(eye(photonic_size), s_z);
sigma_x = kron(eye(photonic_size), s_x);

a_dag_sigma_minus = kron(a_dag, s_0) * kron(eye(photonic_size), s_minus);
sigma_plus_a_std = kron(eye(photonic_size), s_plus) * kron(a_std, s_0);

H0 = 0.5 * delta * sigma_z + 0.5 * omega * (a_dag_sigma_minus + sigma_plus_a_std);
H1 = sigma_x;

% Dissipation
diss = kron(a_std, s_0);

qjx_path = 'D:/Work/os_d/source/cpp/QJX/QJX';
qjx_suffix = sprintf('setup(9_7_1)_rnd(1_1000000)_nph(1)_ampl(%0.4f)_freq(%0.4f)_D(%0.4f)_O(%0.4f)_gamma(%0.4f)_lpn(-1_-2.0000_-2.0000_-2.0000)', ...
    amplitude, ...
    freq, ...
    delta, ...
    omega, ...
    gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn = sprintf('%s/hamiltonian_%s.txt', qjx_path, qjx_suffix);
qjx_data = importdata(fn);
qjx_H = zeros(N, N);
for s_id_1 = 1:N
    for s_id_2 = 1:N
        index = (s_id_1 - 1) * N + s_id_2;
        qjx_H(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
    end
end

H_check = norm(qjx_H - H0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Dissipators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = sprintf('%s/dissipator_%d_%s.txt', qjx_path, 0, qjx_suffix);
qjx_data = importdata(fn);
qjx_diss = zeros(N, N);
for s_id_1 = 1:N
    for s_id_2 = 1:N
        index = (s_id_1 - 1) * N + s_id_2;
        qjx_diss(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
    end
end

diss_check = norm(diss - qjx_diss)
