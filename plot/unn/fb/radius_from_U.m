clear all;

drt = 1;
N = 501;
E = 0;
J = -1;

U_begin = 0.005;
U_shift = 0.005;
U_num = 40;

g = 0.1;
A = -3.4;
omega = 1;
seed = 1;

np = 100;

nu_size = 100;

nu_step = pi / nu_size;

phi = [pi / 2];
phis = phi;
nus = linspace(0, pi, nu_size)';

data_path = '../../../data/cluster/unn';

pks_lim = 0.0001;

warning('off', 'all');

Us = zeros(U_num, 1);

husimis = zeros(U_num, nu_size);

rads = zeros(U_num, 1);

for U_id = 1 : U_num
    
    U = U_begin + (U_id-1) * U_shift
    Us(U_id) = U;
    
    local_path = sprintf('np_%d/drt_%d/N_%d/E0_%0.4f/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_%d', ...
        np, ...
        drt, ...
        N-1, ...
        E, ...
        J, ...
        U, ...
        g, ...
        A, ...
        omega, ...
        seed);
    
    fn = sprintf('%s/%s/rho.txt', data_path, local_path);
    rho_data = importdata(fn);
    
    rho = zeros(N);
    for s_id = 1 : size(rho_data, 1)
        curr_row = rho_data(s_id, 1);
        curr_col = rho_data(s_id, 2);
        rho(curr_row, curr_col) = rho_data(s_id, 3) + sqrt(-1) * rho_data(s_id, 4);
    end
    
    tic
    hus = husimi(nus, phis, rho);
    toc
    
    husimis(U_id, :) = hus;
    
    [pks, locs] = findpeaks(abs(hus));
    
    del_ids = [];
    for i = 1:size(pks, 1)
        if pks(i) < pks_lim
            del_ids = vertcat(del_ids, i);
        end
    end
    pks(del_ids) = [];
    locs(del_ids) = [];
    
    if size(pks, 1) == 1 
        rads(U_id) = 0;
    else
        rads(U_id) = abs(locs(1) - locs(2)) * nu_step;
    end
    
end

fig = figure;
propertyeditor(fig);

hLine = plot(Us, rads);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$R$', 'Interpreter', 'latex');
hold all;

