clear all;

home_figures_path = '/home/yusipov/Work/os_d/figures';

data_path = '/data/biophys/yusipov/os_d';
prefix = 'mf_results';

data_path = sprintf('%s/%s', data_path, prefix);

task = 2;
U = 2.8;
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

epss_log_diff = diff(log10(epss));
cis_log_diff = diff(log10(cis));

local_diff = zeros(eps_num - 1, 1);
for i = 1:eps_num-1
	local_diff(i) = cis_log_diff(i) / epss_log_diff(i);
end

epss_log_scale = log10(epss(1:end-1));

hLine = plot(epss_log_scale, local_diff, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$\log_{10}(\epsilon)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('local diff CD', 'Interpreter', 'latex');


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

savefig(sprintf('%s/mf_cd_ld_%s.fig', home_figures_path, fn_suffix));
