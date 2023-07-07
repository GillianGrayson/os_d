clear all;

N = 2;
N2 = N * N;

delta = 1;

% Term with modulation. Currently only first spin is modulated
amplitude = 0.025;
freq = 0.06;
T  = 2.0 * pi / freq;

% dissipation
gamma = 0.05;

% Matrixes
sigma_0 = eye(2);
sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];
sigma_minus = [0 0; 1 0];

% Hamiltonian
H0 = delta * 0.5 * sigma_z;
H1 = sigma_x;

% Dissipation
disses{1} = sigma_minus;

qjx_path = 'D:/Work/os_d/source/cpp/QJX/QJX';
qjx_suffix = sprintf('setup(8_7_1)_rnd(1_1000000)_n(1)_ampl(%0.4f)_freq(%0.4f)_D(%0.4f)_J(1.0000)_gamma(%0.4f)_lpn(-1_-2.0000_-2.0000_-2.0000)', ...
    amplitude, ...
    freq, ...
    delta, ...
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
for diss_id = 1:1
    fn = sprintf('%s/dissipator_%d_%s.txt', qjx_path, diss_id - 1, qjx_suffix);
    qjx_data = importdata(fn);
    qjx_diss = zeros(N, N);
    for s_id_1 = 1:N
        for s_id_2 = 1:N
            index = (s_id_1 - 1) * N + s_id_2;
            qjx_diss(s_id_1, s_id_2) = qjx_data(index, 1) + 1i * qjx_data(index, 2);
        end
    end
    
    diss_check = norm(disses{diss_id} - qjx_diss)
end
