clear all;

data_path = '../../../source/cpp/MF/MF';

task = 0;
U = 0.12;
seed_start = 0;
seed_num = 1;
path = '';
mt = 1;
num_steps = 1000;
npt = 100;
np = 10000;
E = 0.0;
A = 3.4 ;
omega = 1.0;
phase = 0.0;
gamma = 0.1;
J = 1.0;

nu = [];
phi = [];

for seed = seed_start : seed_start + seed_num-1
    
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
    
    nu_curr = data(2:end, 1);
    phi_curr = data(2:end, 2);
    
    for p_id = 1 : size(nu_curr, 1)
        
        nu_curr(p_id) = mod(nu_curr(p_id), 2 * pi);
        if (nu_curr(p_id) > pi && nu_curr(p_id) < 2.0 * pi)
            nu_curr(p_id) = 2.0 * pi - nu_curr(p_id);
            phi_curr(p_id) = phi_curr(p_id) + pi;
        end
        phi_curr(p_id) = mod(phi_curr(p_id), 2 * pi);
    end
    
    nu = vertcat(nu, nu_curr);
    phi = vertcat(phi, phi_curr);
    
end



fig = figure;
hLine = scatter(nu, phi);
set(gca, 'FontSize', 30);
xlabel('$\nu$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\phi$', 'Interpreter', 'latex');