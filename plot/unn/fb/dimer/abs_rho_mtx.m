clear all;

drt = 1;
N = 501;
E = 0;
J = -1;

U = 0.1125;

g = 0.1;
A = -3.4;
omega = 1;
seed = 1;

np = 100;

data_path = '../../../../data/cluster/unn/fb';

states = linspace(1, N, N)';

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
sprs = importdata(fn);
mtx = zeros(N);

for s_id = 1 : size(sprs, 1)
    curr_row = sprs(s_id, 1);
    curr_col = sprs(s_id, 2);
    mtx(curr_row, curr_col) = sprs(s_id, 3) + sqrt(-1) * sprs(s_id, 4);
end

fig = figure;
propertyeditor(fig);

hLine = imagesc(states, states, abs(mtx));
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$m$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
h.Label.Interpreter = 'latex';
title(h, '$|\rho_{n,m}|$', 'FontSize', 33, 'Interpreter', 'latex');
set(gca,'YDir','normal');
hold all;

fig = figure;
propertyeditor(fig);

diag_mtx = diag(mtx);

hLine = plot(states, diag_mtx);
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\rho_{n,n}$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);

