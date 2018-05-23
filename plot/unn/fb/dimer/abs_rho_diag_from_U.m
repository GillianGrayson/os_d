clear all;

drt = 0;
N = 101;
E = 0;
J = 1;

U_start = 0.005;
U_shift = 0.005;
U_num = 200;

g = 0.1;
A = -1.5;
omega = 1;
seed = 1;

np = 1000;

data_path = '../../../../data/cluster/unn';

Us = zeros(U_num, 1);
Ns = linspace(0, 1, N);

abs_rho_diag = zeros(U_num, N);

for U_id = 1:U_num
    
    U_id = U_id
    
    Us(U_id) = U_start + (U_id-1) * U_shift;
    
   local_path = sprintf('np_%d/drt_%d/N_%d/E0_%0.4f/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_%d', ...
        np, ...
        drt, ...
        N-1, ...
        E, ...
        J, ...
        Us(U_id), ...
        g, ...
        A, ...
        omega, ...
        seed);
    
    fn = sprintf('%s/%s/rho_diag.txt', data_path, local_path);
    abs_rho_diag_curr = importdata(fn);
    
    abs_rho_diag(U_id, :) = abs_rho_diag_curr / max(abs_rho_diag_curr);
    abs_rho_diag(U_id, :) = abs_rho_diag_curr;
    
end

fig = figure;
propertyeditor(fig);

hLine = imagesc(Us, Ns, abs_rho_diag');
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
ylim([0 1])
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
h.Label.Interpreter = 'latex';
title(h, '$\rho_{n,n}$', 'FontSize', 33);
set(gca,'YDir','normal');
hold all;

