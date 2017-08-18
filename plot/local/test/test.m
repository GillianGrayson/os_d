clear all;

N = 100;
delta_max = 0.1;

basic_tr_id = 0;
var_tr_id = 1;

path_to_data = '../../../source/cpp/QJ/QJ';

periods = importdata(sprintf('%s/periods_evo.txt', path_to_data));
mc_basic = importdata(sprintf('%s/mean_evo_trajectory_%d.txt', path_to_data, basic_tr_id));
mc_var = importdata(sprintf('%s/mean_evo_trajectory_%d.txt', path_to_data, var_tr_id));
mc_var_real = importdata(sprintf('%s/mean_real_evo_trajectory_%d.txt', path_to_data, var_tr_id));
lam_var = importdata(sprintf('%s/lambda_evo_trajectory_%d.txt', path_to_data, var_tr_id));

t_jump_basic = importdata(sprintf('%s/jump_times_trajectory_%d.txt', path_to_data, basic_tr_id));
t_jump_var = importdata(sprintf('%s/jump_times_trajectory_%d.txt', path_to_data, var_tr_id));

mc_jump_basic = importdata(sprintf('%s/jump_means_trajectory_%d.txt', path_to_data, basic_tr_id));
mc_jump_var = importdata(sprintf('%s/jump_means_trajectory_%d.txt', path_to_data, var_tr_id));



abs_diff = abs(mc_basic - mc_var);

num_dumps = size(periods, 1);

renorms = zeros(num_dumps, 1);
for dump_id = 1:num_dumps
   curr_delta = abs(mc_var_real(dump_id) - mc_basic(dump_id)) / N;
   if  ((curr_delta > delta_max) || (curr_delta < 1.0e-13))
       renorms(dump_id) = 1;
   end
end

figure
hLine = plot(periods, mc_basic,  'LineWidth', 2);
legend(hLine, {'$n_b$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;
hLine = plot(periods, mc_var,  'LineWidth', 2);
legend(hLine, {'$n_v$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;
hLine = plot(periods, mc_var_real,'.-', 'LineWidth', 0.5);
legend(hLine, {'$n_v$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;

figure
hLine = plot(periods, lam_var,  'LineWidth', 2);
legend(hLine, {'$\lambda$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;

hLine = plot(periods, renorms, 'o', 'LineStyle', 'none',  'LineWidth', 1, 'MarkerEdgeColor', [0 0 0]);
legend(hLine, {'$vars$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;

figure
hLine = plot(t_jump_basic / (2*pi), mc_jump_basic, 'LineWidth', 2);
legend(hLine, {'$basic$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;

hLine = plot(t_jump_var / (2*pi), mc_jump_var, 'LineWidth', 2);
legend(hLine, {'$var$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;

hLine = plot(periods - 1, mc_basic, 'o', 'LineStyle', 'none',  'LineWidth', 1);
legend(hLine, {'$basic$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;

hLine = plot(periods - 1, mc_var_real, 'o', 'LineStyle', 'none',  'LineWidth', 1);
legend(hLine, {'$var$'}, 'Interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t, periods$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
hold all;


