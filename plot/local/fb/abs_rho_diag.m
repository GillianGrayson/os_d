clear all;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/rho_diag.txt', data_path);
rho_diag_data = importdata(fn);

figure;
plot(rho_diag_data, 'LineWidth', 2);set(gca, 'FontSize', 30);
xlabel('$n$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\rho_{n,n}$', 'Interpreter', 'latex');
hold all;

