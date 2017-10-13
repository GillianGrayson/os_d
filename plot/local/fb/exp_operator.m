clear all;

N = 5;

NS = (N+1) * (N+1) - 1;

T = 2 * pi;

fig = figure;
propertyeditor(fig);

circle(0, 0, 1)
hold all;
xlim([-2 2])
ylim([-2 2])

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

exp_prop_small = expm(Gs * T);
evals_small = eig(exp_prop_small);

scatter(real(evals_small), imag(evals_small), '.');
hold all

fn = sprintf('%s/Ks.txt', data_path);
data = 0;
data = importdata(fn);
d_size = size(data, 1);

Ks = zeros(NS);

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

exp_prop_big = expm(op_big * T);
evals_big = eig(exp_prop_big);

scatter(real(evals_big), imag(evals_big), 'LineWidth', 2);
hold all


function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end

