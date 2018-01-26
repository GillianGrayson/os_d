clear all;

data_path = '../../../source/cpp/MF/MF';

task = 2;
mt = 0;
omega = 1.0;
phase = 0.0;
gamma = 0.1;
J = 1.0;
E = 1.0;
A = 1.5 ;

U = 0.1;

seed = 0;

fn_suffix = sprintf('t(%d)_mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed(%d).txt', ...
    task, ...
    mt, ...
    omega, ...
    phase, ...
    gamma, ...
    J, ...
    E, ...
    A, ...
    U, ...
    seed);

fn = sprintf('%s/data_%s', data_path, fn_suffix);
data = importdata(fn);

y = data(:,1);
y = (cos(y) + 1) * 0.5;

x = data(:,end) / (2*pi);

fig = figure;
hLine = plot(x, y);
title_str = sprintf('$E=%0.2f$ $U=%0.2f$ $J=%0.2f$ ', E, U, J);
title(title_str, 'FontSize', 33, 'interpreter','latex');
set(gca, 'FontSize', 30);
xlabel('$t/T$', 'Interpreter', 'latex');
xlim([x(1) x(end)])
set(gca, 'FontSize', 30);
ylabel('$y$', 'Interpreter', 'latex');
hold all;

propertyeditor(fig)
