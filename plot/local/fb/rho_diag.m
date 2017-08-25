clear all;

N = 1000;
J =  -1; 
E0 = 0; 
U = 0.5; 
g = 0.1; 
A0 = -3.4; 
omega = 1;

data_path = '../../data/cluster/unn';

fn_path = sprintf('N_%d/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_1', ...
    N, ...
    J, ...
    U, ...
    g, ...
    A0, ...
    omega);
fn = sprintf('%s/%s/rho_diag.txt', data_path, fn_path);
rd = importdata(fn);

Ns = linspace(1, N+1, N+1);

figure;
hLine = plot(Ns, rd, 'LineWidth', 2);
title(sprintf('J=%0.1f E_0=%0.1f A=%0.1f \\omega=%0.1f', J, E0, A0, omega));
set(gca, 'FontSize', 30);
xlabel('$N$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$|\rho_{n,n}|$', 'Interpreter', 'latex');
hold all;