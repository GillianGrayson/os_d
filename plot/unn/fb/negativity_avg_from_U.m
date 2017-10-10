clear all;

drt = 1;
N = 101;
E = -1;
J = -1;

U_start = 0.005;
U_shift = 0.005;
U_num = 150;

g = 0.1;
A = -3.4;
omega = 1;
seed = 1;

data_path = '../../../data/cluster/unn';

Us = zeros(U_num, 1);
neg_avg = zeros(U_num, 1);

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
    
    fn = sprintf('%s/%s/negativity_avg.txt', data_path, local_path);
    neg_avg_curr = importdata(fn);
    
    neg_avg(U_id) = neg_avg_curr(1);
    
end

fig = figure;
propertyeditor(fig);

hLine = plot(Us, neg_avg, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\bar N$', 'Interpreter', 'latex');
hold all;

