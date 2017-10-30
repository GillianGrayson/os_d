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

phi_size = 500;
nu_size = 500;

phis = linspace(0, 2*pi, phi_size);
nus = linspace(0, pi, nu_size);

data_path = '../../../data/cluster/unn';


local_path = sprintf('drt_%d/N_%d/E0_%0.4f/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_%d', ...
    drt, ...
    N-1, ...
    E, ...
    J, ...
    Us(U_id), ...
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

hus = husimi(nus, phis, rho);
    

fig = figure;
propertyeditor(fig);

hLine = imagesc(nus, phis, hus);
set(gca, 'FontSize', 30);
xlabel('$\nu$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\phi$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, 'H', 'FontSize', 33);
set(gca,'YDir','normal');
ylim([Ns(1) - 0.5 Ns(end) + 0.5]);
hold all;

