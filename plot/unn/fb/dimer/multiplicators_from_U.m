clear all;

drt = 1;
N = 26;
E = -1;
J = -1;

U_start = 0.005;
U_shift = 0.005;
U_num = 150;

g = 0.1;
A = -1.5;
omega = 1;
seed = 1;

mults_num = N*N;
mult_id = 1;

data_path = '../../../data/cluster/unn';

Us = zeros(U_num, 1);
abs_mults = zeros(U_num, 1);

for U_id = 1:U_num
    
    Us(U_id) = U_start + (U_id-1) * U_shift;
    
    local_path = sprintf('drt_%d/N_%d/E0_%0.4f/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_%d', ...
        drt, ...
        N-1, ...
        E, ...
        J, ...
        Us(U_id), ...
        g, ...
        A, ...
        omega, ...
        seed);
    
    fn = sprintf('%s/%s/multiplicators.txt', data_path, local_path);
    curr_data = importdata(fn);
    
    curr_mults = curr_data(:, 1) + sqrt(-1) * curr_data(:, 2);
    curr_abs_mults = abs(curr_mults);
    
    curr_abs_mults = sort(curr_abs_mults, 'descend');
    
    abs_mults(U_id) = curr_abs_mults(mult_id);
end

fig = figure;
propertyeditor(fig);

hLine = plot(Us, abs_mults, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\mu$', 'Interpreter', 'latex');
hold all;

