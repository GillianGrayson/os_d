clear all;

drt = 1;
N = 501;
E = 0;
J = -1;

U = 0.15;

g = 0.1;
A = -3.4;
omega = 1;
seed = 1;

np = 100;

phi_size = 100;
nu_size = 100;

phis = linspace(0, 2*pi, phi_size)';
nus = linspace(0, pi, nu_size)';

data_path = '../../../data/cluster/unn';

warning('off', 'all');

local_path = sprintf('np_%d/drt_%d/N_%d/E0_%0.4f/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_%d', ...
    np, ...
    drt, ...
    N-1, ...
    E, ...
    J, ...
    U, ...
    g, ...
    A, ...
    omega, ...
    seed);

fn = sprintf('%s/%s/rho.txt', data_path, local_path);
rho_data = importdata(fn);

rho = zeros(N);
for s_id = 1 : size(rho_data, 1)
    curr_row = rho_data(s_id, 1);
    curr_col = rho_data(s_id, 2);
    rho(curr_row, curr_col) = rho_data(s_id, 3) + sqrt(-1) * rho_data(s_id, 4);
end

tic
hus = husimi(nus, phis, rho);
toc 

fig = figure;
propertyeditor(fig);

hLine = imagesc(nus, phis, real(hus'));
set(gca, 'FontSize', 30);
xlabel('$\theta$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\phi$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'H', 'FontSize', 33);
set(gca,'YDir','normal');
hold all;

