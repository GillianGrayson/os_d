clear all;

basic_tr_id = 0;
var_tr_id = 1;

periods = importdata(sprintf('periods_evo.txt'));
mc_basic = importdata(sprintf('mean_evo_trajectory_%d.txt', basic_tr_id));
mc_var = importdata(sprintf('mean_evo_trajectory_%d.txt', var_tr_id));
lam_var = importdata(sprintf('lambda_evo_trajectory_%d.txt', var_tr_id));

abs_diff = abs(mc_basic - mc_var);

figure
hLine = plot(periods, mc_basic,  'LineWidth', 2);
legend({'$n_b$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;

figure
hLine = plot(periods, abs_diff,  'LineWidth', 2);
legend({'$|n_b - n_v|$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;

hLine = plot(periods, lam_var,  'LineWidth', 2);
legend({'$\lambda$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;


