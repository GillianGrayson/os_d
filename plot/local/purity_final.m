clear all;
 
J =  -1; 
E0 = 0; 
U = 0.5; 
g = 0.1; 
A0 = -3.4; 
omega = 1;

data_path = '../../data/cluster/unn';

N_start = 100;
N_num = 10;
N_shift = 100;

purities_final = zeros(N_num, 1);
purities_final_avg = zeros(N_num, 1);
Ns = zeros(N_num, 1);

for N_id = 1:N_num
    
    N = N_start + (N_id - 1) * N_shift;
    Ns(N_id) = N;
    
    fn_path = sprintf('N_%d/J_%0.4f/U_%0.4f/g_%0.4f/A0_%0.4f/omega_%0.4f/seed_1', ...
        N, ...
        J, ...
        U, ...
        g, ...
        A0, ...
        omega);
    fn = sprintf('%s/%s/purity_final.txt', data_path, fn_path);
    pur_fin_data = importdata(fn);
    pur_fin = pur_fin_data(1);
    purities_final(N_id) = pur_fin;

    fn = sprintf('%s/%s/purity_avg.txt', data_path, fn_path);
    pur_avg_data = importdata(fn);
    pur_avg = pur_avg_data(1);
    purities_final_avg(N_id) = pur_avg;
end

hLine = plot(Ns, purities_final, 'LineWidth', 2);
title(sprintf('J=%0.1f E_0=%0.1f A=%0.1f \\omega=%0.1f', J, E0, A0, omega));
legend(hLine, sprintf('last time'));
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$purity$', 'Interpreter', 'latex');
hold all;

hLine = plot(Ns, purities_final_avg, 'LineWidth', 2);
legend(hLine, sprintf('avg last period'));
set(gca, 'FontSize', 30);
xlabel('$t$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$purity$', 'Interpreter', 'latex');
hold all;