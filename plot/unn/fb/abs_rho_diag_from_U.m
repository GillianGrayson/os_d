clear all;


E:\Work\os_d\data\cluster\unn\drt_1\N_100\E0_-1.0000\J_-1.0000\U_0.0300\g_0.1000\A0_-3.4000\omega_1.0000\seed_1

drt = 1;
N = 100;
E = -1;
J = -1;

U_start = 0.03;
U_shift = 0.03;
U_num = 100;

g = 0.1;
A = -3.4;
omega = 1;
seed = 1;

data_path = '../../../data/cluster/unn';

Us = zeros(U_num, 1);
Ns = linspace(1, N, N);

abs_rho_diag = zeros(U_num, N);

for U_id = 1:U_num
    
    Us(U_id) = U_start + (U_id-1) * U_shift;
    
    local_path = sprintf('drt_%d/N_%d/E0_%0.4f/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_%d', ...
        drt, ...
        N, ...
        E, ...
        J, ...
        Us(U_id), ...
        g, ...
        A, ...
        omega, ...
        seed);
    
    fn = sprintf('%s/%s/rho_diag.txt', data_path, local_path);
    abs_rho_diag_curr = importdata(fn);
    
    abs_rho_diag(U_id, :) = abs_rho_diag_curr;
    
end

figure;
hLine = imagesc(Us, Ns, abs_rho_diag);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\rho_{n,n}', 'FontSize', 33, 'Interpreter', 'latex');
set(gca,'YDir','normal');
ylim([Ns(1) - 0.5 Ns(end) + 0.5]);
hold all;

