clear all;

%fig = open("D:/YandexDisk/Work/os_d/drafts/floquet/two_spins/norms/norm2d_1_ampl(0.2000_0.2000_500)_freq(0.0200_0.0200_500)_phase(0.0000_0.0000_0).fig");
fig = open("D:/YandexDisk/Work/os_d/drafts/floquet/one_spin/norms/mu1_setup(1)_method(1.00e+06_1.00e-04_1.00e-12)_model(1.0000_0.0500)_modulation(1000_1.0000_0.0025_0.0050_500_0.0060_0.0120_500_0.0000).fig");
a = get(gca,'Children');
x = get(a, 'XData');
y = get(a, 'YData');
z = get(a, 'CData');
z(isinf(z)|isnan(z)) = -12;
close(fig)

%fig = open("D:/YandexDisk/Work/os_d/drafts/floquet/two_spins/lambdas/lambdas_traj(50)_tp(100)_obs(100)_ampl(1.0000_1.0000_100)_freq(0.1000_0.1000_100)_D1(1.0000)_D2(1.0000)_J(1.0000)_gamma(0.0500)_lpn(-1_-1.0000_-1.0000_-1.0000).fig");
fig = open("D:/YandexDisk/Work/os_d/drafts/floquet/many_spins/lambdas/lambdas_traj(50)_tp(100)_obs(100)_n(1)_ampl(0.0250_0.0500_50)_freq(0.0600_0.1200_50)_D(1.0000)_J(1.0000)_gamma(0.0500)_lpn(-1_-1.0000_-1.0000_-1.0000).fig");
hold all;
level = [-11, -10]; 
[C, h] = contour(x, y, z, level, 'LineWidth', 2);
w = h.EdgeColor;
h.EdgeColor = 'k';