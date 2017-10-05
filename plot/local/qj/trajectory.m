clear all;
imag1 = sqrt(-1);

trajectory_id = 0;

N = 100;
E = 0.0;
A = 0.0;
U = 1.0;

ch_name = 'energy';

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

file_name = sprintf('%s/%s_evo_trajectory_%d.txt', qj_data_path, ch_name, trajectory_id);
ch = importdata(file_name);


fig_tr = figure;


subplot(2,1,1);
hLine = imagesc(periods, states, abs_diag_rho_evol);
title_str = sprintf('$N=%d$ $E=%0.2f$ $A=%0.2f$ $\\mathbf{U=%0.2f}$', N, E, A, U);
title(title_str, 'FontSize', 33, 'interpreter','latex');
set(gca, 'FontSize', 30);
%xlabel('$t$', 'Interpreter', 'latex');
set(gca, 'XAxisLocation', 'top');
set(gca,'xticklabel',[]);
xlim([periods(1) periods(end)])
set(gca, 'FontSize', 30);
ylabel('$n$', 'Interpreter', 'latex');
colormap hot;
h = colorbar;
set(gca, 'FontSize', 30);
title(h, '$\rho_{n,n}$', 'FontSize', 33, 'interpreter','latex');
set(gca,'YDir','normal');
set(gca, 'Position', [0.15 0.50 0.70 0.40])
hold all;

subplot(2,1,2);
plot(periods, ch, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
xlim([periods(1) periods(end)])
set(gca, 'FontSize', 30);
ylabel(sprintf('%s', ch_name), 'Interpreter', 'latex');
set(gca, 'Position', [0.15 0.15 0.70 0.35])

propertyeditor(fig_tr)


