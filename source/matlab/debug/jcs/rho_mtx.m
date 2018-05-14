clear all;

cpp_path = '../../../cpp/os_jcs/os_jcs';

N = 50; 

fn = sprintf('%s/rho.txt', cpp_path);
data = importdata(fn);
rho = zeros(N);
for s_id = 1 : size(data, 1)
    curr_row = data(s_id, 1);
    curr_col = data(s_id, 2);
    rho(curr_row, curr_col) = data(s_id, 3) + sqrt(-1) * data(s_id, 4);
end
figure
plot(abs(diag(rho)))