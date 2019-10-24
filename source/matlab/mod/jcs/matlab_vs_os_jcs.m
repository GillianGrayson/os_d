clear all;

cpp_path = '../../../cpp/os_jcs/os_jcs';

alfa=5;

N = 50; 
a_std = zeros(N,N);
for n=1:(N-1)
    a_std(n,n+1)=sqrt(n);
end
a_dag=a_std';

%??????????? ???????????? H(t)=H1+F(t)*H2
H_0=1/2*(1/(alfa*alfa*alfa))*a_dag*a_dag*a_std*a_std;
H_1=i*(a_dag-a_std);

fn = sprintf('%s/H0.txt', cpp_path);
data = importdata(fn);
H_0_cpp = zeros(N);
for s_id = 1 : size(data, 1)
    curr_row = data(s_id, 1);
    curr_col = data(s_id, 2);
    H_0_cpp(curr_row, curr_col) = data(s_id, 3) + sqrt(-1) * data(s_id, 4);
end

H_0_diff = max(max(H_0_cpp - H_0))

fn = sprintf('%s/H1.txt', cpp_path);
data = importdata(fn);
H_1_cpp = zeros(N);
for s_id = 1 : size(data, 1)
    curr_row = data(s_id, 1);
    curr_col = data(s_id, 2);
    H_1_cpp(curr_row, curr_col) = data(s_id, 3) + sqrt(-1) * data(s_id, 4);
end

H_1_diff = max(max(H_1_cpp - H_1))
