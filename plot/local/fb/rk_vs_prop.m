clear all;

N = 5;

NS = (N+1) * (N+1) - 1;

T = 2 * pi;
NT = 11;


data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/Gs.txt', data_path);
data = 0;
data = importdata(fn);
d_size = size(data, 1);


Gs = zeros(NS, NS);
for d_id = 1:d_size
    curr_row = data(d_id, 1);
    curr_col = data(d_id, 2);
    Gs(curr_row, curr_col) = data(d_id, 3) + sqrt(-1) * data(d_id, 4);
end


fn = sprintf('%s/Ks.txt', data_path);
data = 0;
data = importdata(fn);
d_size = size(data, 1);

Ks = zeros(NS, 1);

for d_id = 1:d_size
    Ks(d_id) = data(d_id, 1) + sqrt(-1) * data(d_id, 2);
end

op_big = zeros(NS + 1, NS + 1);

for s_id = 1:NS
    op_big(s_id + 1, 1) = -Ks(s_id);
end

for s_id_1 = 1 : NS
    for s_id_2 = 1 : NS
        op_big(s_id_1 + 1, s_id_2 + 1) = Gs(s_id_1, s_id_2);
    end
end

exp_prop_big = expm(op_big * NT * T);
exp_prop_big(1, :) = [];
exp_prop_big(:, 1) = [];

fn = sprintf('%s/rho_f_init.txt', data_path);
data = 0;
data = importdata(fn);
rho_f_init = data(:, 1);

fn = sprintf('%s/rho_f_fin.txt', data_path);
data = 0;
data = importdata(fn);
rho_f_fin = data(:, 1);
