clear all;

N = 21;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/rho.txt', data_path);
rho_data = importdata(fn);

d_size = size(rho_data, 1);

states = linspace(1, N, N);

curr_rho = zeros(N, N);
for d_id = 1:d_size
    curr_row = rho_data(d_id, 1);
    curr_col = rho_data(d_id, 2);
    curr_rho(curr_row, curr_col) = rho_data(d_id, 3) + sqrt(-1) * rho_data(d_id, 4);
end

fig = figure;
hLine = imagesc(states, states, abs(curr_rho));
set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$m$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$|\rho_{n,m}|$', 'Interpreter', 'latex');
set(gca,'YDir','normal');

rho_2 = curr_rho * curr_rho;
tr = trace(rho_2)

propertyeditor('on')

