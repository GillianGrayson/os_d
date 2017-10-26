clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'mf_results';

data_path = sprintf('%s/%s', data_path, prefix);

task = 2;
U = 0.3;
seed_start = 1;
seed_num = 1;
mt = 0;
num_steps = 10000;
npt = 2000;
np = 100000;
E = 1.0;
A = 1.5 ;
omega = 1.0;
phase = 0.0;
gamma = 0.1;
J = 1.0;

eps_start = 1.0e-8;
eps_per_decade = 10;
eps_num = 81;

begin_fit_id = 31;
end_fit_id = 46;

fig = figure;

m = 4;

epss = zeros(eps_num, 1);
cis = zeros(eps_num, 1);

for eps_id = 1:eps_num
    
    eps_id = eps_id
    
    eps = eps_start * 10.0^((eps_id-1)/eps_per_decade);
    
    epss(eps_id) = eps;
    
    for seed = 1:seed_num
        
        path_to_folder = sprintf('%s/task_%d/mt_%d/omega_%0.4f/phase_%0.4f/g_%0.4f/J_%0.4f/E_%0.4f/A_%0.4f/U_%0.4f/eps_%0.10f/m_%d/seed_%d', ...
            data_path, ...
            task, ...
            mt, ...
            omega, ...
            phase, ...
            gamma, ...
            J, ...
            E, ...
            A, ...
            U, ...
            eps, ...
            m, ...
            seed);
        
        fn_suffix = sprintf('eps(%0.10f)_m(%d)_mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed(%d).txt', ...
            eps, ...
            m, ...
            mt, ...
            omega, ...
            phase, ...
            gamma, ...
            J, ...
            E, ...
            A, ...
            U, ...
            seed);
        
        fn = sprintf('%s/ci_%s', path_to_folder, fn_suffix);
        data = importdata(fn);
        
        cis(eps_id) = cis(eps_id) + data;
    end
    
    cis(eps_id) = cis(eps_id) / seed_num;
    
end

epss_log = log10(epss);
cis_log = log10(cis);

fit_epss_log = epss_log(begin_fit_id:end_fit_id);
fit_cis_log = cis_log(begin_fit_id:end_fit_id);

[lin_fit, gof] = fit(fit_epss_log, fit_cis_log, 'poly1');
lin_formula = formula(lin_fit);
coeffs_lin_fit = coeffvalues(lin_fit);
D = coeffs_lin_fit(1);
conf_int_coeffs = confint(lin_fit, 0.95);
conf_int_D = conf_int_coeffs(:,1)

fitted_cis_log(:) = D * fit_epss_log(:) + coeffs_lin_fit(2);

hLine = plot(epss_log, cis_log, 'LineWidth', 2);
title(sprintf('D=%0.4f left=%0.4f right=%0.4f', D, conf_int_D(1), conf_int_D(2)));
legend(hLine, sprintf('m=%d', m));
set(gca, 'FontSize', 30);
xlabel('$\log_{10}(\epsilon)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\log_{10}(C)$', 'Interpreter', 'latex');
hold all;
hLine = plot(fit_epss_log, fitted_cis_log, 'LineWidth', 2);
legend(hLine, sprintf('m=%d fitted', m));


fn_suffix = sprintf('mt(%d)_omega(%0.4f)_phase(%0.4f)_g(%0.4f)_J(%0.4f)_E(%0.4f)_A(%0.4f)_U(%0.4f)_seed_num(%d)', ...
    mt, ...
    omega, ...
    phase, ...
    gamma, ...
    J, ...
    E, ...
    A, ...
    U, ...
    seed_num);

savefig(sprintf('%s/mf_cd_fitted_%s.fig', home_figures_path, fn_suffix));
