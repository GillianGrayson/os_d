clear all;

N = 5;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/rho_34_init.txt', data_path);
rho_data = importdata(fn);

d_size = size(rho_data, 1);

states = linspace(1, N, N);

mtx = zeros(N + 1, N + 1);
for d_id = 1:d_size
    curr_row = rho_data(d_id, 1);
    curr_col = rho_data(d_id, 2);
    mtx(curr_row, curr_col) = rho_data(d_id, 3) + sqrt(-1) * rho_data(d_id, 4);
end

mtx = mtx - eye(N + 1) / (N+1);

fig = figure;
hLine = imagesc(states, states, mtx);
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$m$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$|\rho_{n,m}|$', 'Interpreter', 'latex');
set(gca,'YDir','normal');

tr = trace(mtx)

propertyeditor('on')

