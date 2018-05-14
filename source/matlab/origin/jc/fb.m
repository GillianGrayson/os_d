clear all;
tic

imag1=sqrt(-1);

data_path = '../../../cpp/os_jcs/os_jcs';

N = 9;
M = N+1;
prm_alpha = 5;
drv_ampl = 3.2;
g = 0.1;
gamma = 0.25 / prm_alpha;

a_std = zeros(M,M);
for n = 1:(M-1)
    a_std(n,n+1) = sqrt(n);
end
a_dag = a_std';

V = a_std;

H1 = 1/2 * (1/(prm_alpha * prm_alpha * prm_alpha)) * a_dag*a_dag*a_std*a_std;
H_drv = i*(a_dag - a_std);
H0 = H1 + drv_ampl * H_drv;

num_steps_t_0 = 1000;
num_steps_t_1 = 1000;
num_periods_trans = 100;
num_periods_obser = 1;

t_0 = 0.98 * prm_alpha;
t_1 = 1.0 * prm_alpha;

%%% f-basis
k=1;
F{1}=sparse(eye(N+1));
for i=1:N+1
    for j=i+1:N+1
        k=k+1;
        F{k}=sparse([i j],[j i],[1 1]/sqrt(2),N+1,N+1);
        k=k+1;
        F{k}=sparse([i j],[j i],-imag1*[1 -1]/sqrt(2),N+1,N+1);
    end
end

for i=1:N
    k=k+1;
    temp(1:i)=ones(1,i);
    temp(i+1)=-i;
    F{k}=sparse([1:i+1],[1:i+1],temp/sqrt(i*(i+1)),N+1,N+1);
end

H0s = sparse(H0);
H_drv_s = sparse(H0);
H1s = sparse(H1);

%% H0s Checking
fn = sprintf('%s/H_0.txt', ...
    data_path);
sprs = importdata(fn);

H0_cpp = zeros(size(H_drv_s, 1));

for s_id = 1 : size(sprs, 1)
    curr_row = sprs(s_id, 1);
    curr_col = sprs(s_id, 2);
    H0_cpp(curr_row, curr_col) = sprs(s_id, 3) + sqrt(-1) * sprs(s_id, 4);
end

H0_diff = max(max(abs(H0_cpp - H1)))

%% H1s Checking
fn = sprintf('%s/H_1.txt', ...
    data_path);
sprs = importdata(fn);

H1_cpp = zeros(size(H_drv_s, 1));

for s_id = 1 : size(sprs, 1)
    curr_row = sprs(s_id, 1);
    curr_col = sprs(s_id, 2);
    H1_cpp(curr_row, curr_col) = sprs(s_id, 3) + sqrt(-1) * sprs(s_id, 4);
end

H1_diff = max(max(abs(H1_cpp - H_drv)))


for i=1:(N+1)^2-1
    
    h0s(i)=trace(H0s*F{i+1});
    h1s(i)=trace(H1s*F{i+1});
    
end

a_std=zeros((N+1)^2-1);

A1=V;
A2=zeros(M);

A1s=sparse(V);
A1=0;
A2s=sparse(A2);
A2=0;

for i=1:(N+1)^2-1
    a1(i)=trace(A1s*F{i+1});
    a2(i)=trace(A2s*F{i+1});
end

for i=1:(N+1)^2-1
    for j=1:(N+1)^2-1
        a_std(i,j)=a_std(i,j)+(a1(i)-imag1*a2(i))*conj(a1(j)-imag1*a2(j));
    end
end


as=sparse(a_std);
a_std=0;

for i=1:(N+1)^2-1
    for j=1:(N+1)^2-1
        for k=1:(N+1)^2-1
            f_temp(j,k)=-imag1*trace(F{k+1}*(F{i+1}*F{j+1}-F{j+1}*F{i+1}));
            d_temp(j,k)=trace(F{k+1}*(F{i+1}*F{j+1}+F{j+1}*F{i+1}));
        end
    end
    fs{i}=sparse(f_temp);
    nzf(i)=nzmax(fs{i});
    ds{i}=sparse(d_temp);
    nzd(i)=nzmax(ds{i});
end
f_temp=0;
d_temp=0;

Q0s=sparse(zeros((N+1)^2-1));
Q1s=sparse(zeros((N+1)^2-1));


for i=1:(N+1)^2-1
    
    Q0s=Q0s+h0s(i)*transpose(fs{i});
    Q1s=Q1s+h1s(i)*transpose(fs{i});
    
end

K=zeros((N+1)^2-1,1);
for j=1:(N+1)^2-1
    k_temp=0;
    for i=1:(N+1)^2-1
        k_temp=k_temp+as(i,:)*fs{i}(:,j);
    end
    K(j)=k_temp;
end

K=imag1/(N+1)*K;

Ks=sparse(K);
Ks=gamma*Ks;
K=0;

R=zeros((N+1)^2-1);
for i=1:(N+1)^2-1
    for k=1:(N+1)^2-1
        R=R+as(i,k)*(fs{k}.'*(fs{i}+imag1*ds{i})+fs{i}.'*conj(fs{k}+imag1*ds{k}));
    end
end

R=-R/4;

Rs=sparse(R);
Rs=gamma*Rs;
R=0;

G0s = Q0s + Rs;
G1s = Q1s + Rs;


if(max(max(abs(imag(G0s))))<1e-15)
    G0s=real(G0s);
end

if(max(max(abs(imag(G1s))))<1e-15)
    G1s=real(G1s);
end

if(max(abs(imag(Ks)))<1e-15)
    Ks=real(Ks);
end

toc

%% G0s Checking
fn = sprintf('%s/G_0_s.txt', ...
    data_path);
sprs = importdata(fn);

G0s_cpp = zeros(size(G0s, 1));

for s_id = 1 : size(sprs, 1)
    curr_row = sprs(s_id, 1);
    curr_col = sprs(s_id, 2);
    G0s_cpp(curr_row, curr_col) = sprs(s_id, 3) + sqrt(-1) * sprs(s_id, 4);
end

G0s_diff = max(max(abs(G0s_cpp - G0s)))

%% G1s Checking
fn = sprintf('%s/G_1_s.txt', ...
    data_path);
sprs = importdata(fn);

G1s_cpp = zeros(size(G0s, 1));

for s_id = 1 : size(sprs, 1)
    curr_row = sprs(s_id, 1);
    curr_col = sprs(s_id, 2);
    G1s_cpp(curr_row, curr_col) = sprs(s_id, 3) + sqrt(-1) * sprs(s_id, 4);
end

G1s_diff = max(max(abs(G1s_cpp - G1s)))


%% Ks Checking
fn = sprintf('%s/Ks_%s.txt', ...
    data_path);
sprs = importdata(fn);
Ks_cpp = -sprs(:, 1);

Ks_diff = max(abs(Ks_cpp - Ks))

RhoF = zeros((N+1)^2-1, 1);

for per_id = 1:num_periods
    
    [t,Y] = ode45('rp_H0', [0 t_0], RhoF);
    RhoF=Y(length(t),:)';
    
    [t,Y] = ode45('rp_H1', [0 t_1], RhoF);
    RhoF=Y(length(t),:)';
end