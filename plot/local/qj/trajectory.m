clear all;
imag1 = sqrt(-1);

trajectory_id = 0;

N = 50;

qj_data_path = '../../../source/cpp/QJ/QJ';


file_name = sprintf('%s/periods_evo.txt', qj_data_path);
periods = importdata(file_name);
periods = periods - 1;
num_dumps = size(periods, 1);

abs_diag_rho_evol = zeros(N, num_dumps);

file_name = sprintf('%s/abs_diag_rho_trajectory_%d.txt', qj_data_path, trajectory_id);
abs_diag_rho_evol_data = importdata(file_name);

for dump_id = 1:num_dumps
    for i = 1:N
        abs_diag_rho_evol(i, dump_id) = abs_diag_rho_evol_data((dump_id-1)*N + i);
    end
end

states = zeros(N,1);
for i = 1:N
    states(i) = i;
end

hLine = imagesc(periods, states, abs_diag_rho_evol);
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '\rho_{n,n}', 'FontSize', 33);
set(gca,'YDir','normal');
hold all;

