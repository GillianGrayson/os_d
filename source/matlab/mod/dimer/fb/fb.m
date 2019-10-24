clear all;

global N J E0 A0 w U g QEs QUs QJs Rs Ks phi0 A Ps Pd

path_to_debug = 'E:/Work/os_d/source/cpp/CQdiss_os_dimer/CQdiss_fbasis';

N=10; % number of particles, system size = N+1
imag1=sqrt(-1);

J=-1; % hopping constant
E0=0.; % bias
A0=-3.4; %driving amplitude
w=1; % driving frequency
phi0=pi*0; % initial phase
U=0.5; % on-site interaction
g=0.1/N; % gamma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating H in F-basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HE=zeros(N+1);
HJ=zeros(N+1);
HU=zeros(N+1);

HE=diag(2*[0:N]-N); % diagonal
HU=diag((2*[0:N]-N).^2); % diagonal
for i=1:N
    HJ(i+1,i)=sqrt((N-i+1)*(i));
    HJ(i,i+1)=sqrt((i)*(N-i+1));
end

HE=HE-trace(HE)/(N+1)*eye(N+1);
HU=HU-trace(HU)/(N+1)*eye(N+1);
HJ=HJ-trace(HJ)/(N+1)*eye(N+1);

H=(E0*HE+U/N*HU+J*HJ);

fn = sprintf('%s/H_base.txt', path_to_debug);
mtx_data = importdata(fn);
mtx = zeros(N+1);
for st = 1:size(mtx_data, 1)
    mtx(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
end
H_diff = max(max((abs(H - mtx))))

HEs=sparse(HE); 
HUs=sparse(HU); 
HJs=sparse(HJ); 

for i=1:(N+1)^2-1
    hE(i)=trace(HEs*F{i+1});
    hU(i)=trace(HUs*F{i+1});
    hJ(i)=trace(HJs*F{i+1});
end

A1s=sparse(diag(2*[0:N]-N));
A1 = full(A1s);
A2=zeros(N+1);
for i=1:N
    A2(i+1,i)=-sqrt((N-i+1)*(i));
    A2(i,i+1)=sqrt((i)*(N-i+1));
end
A2=imag1*A2;
A2s=sparse(A2);

fn = sprintf('%s/A1.txt', path_to_debug);
mtx_data = importdata(fn);
mtx = zeros(N+1);
for st = 1:size(mtx_data, 1)
    mtx(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
end
A1_diff = max(max((abs(full(A1s) - mtx))))

fn = sprintf('%s/A2.txt', path_to_debug);
mtx_data = importdata(fn);
mtx = zeros(N+1);
for st = 1:size(mtx_data, 1)
    mtx(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
end
A2_diff = max(max((abs(full(A2s) - mtx))))

for i=1:(N+1)^2-1
    a1(i)=trace(A1s*F{i+1});
    a2(i)=trace(A2s*F{i+1});
end

fn = sprintf('%s/A1F.txt', path_to_debug);
mtx_data = importdata(fn);
vec = zeros((N+1)^2-1, 1);
for st = 1:size(mtx_data, 1)
    vec(mtx_data(st, 1)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
end
A1F_diff = max(abs(a1' - vec))

fn = sprintf('%s/A2F.txt', path_to_debug);
mtx_data = importdata(fn);
vec = zeros((N+1)^2-1, 1);
for st = 1:size(mtx_data, 1)
    vec(mtx_data(st, 1)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
end
A2F_diff = max(abs(a2' - vec))

for i=1:(N+1)^2-1
    for j=1:(N+1)^2-1
        a(i,j)=(a1(i)-imag1*a2(i))*conj(a1(j)-imag1*a2(j));
    end
end

fn = sprintf('%s/diss_mtx.txt', path_to_debug);
mtx_data = importdata(fn);
mtx = zeros((N+1)^2-1);
for st = 1:size(mtx_data, 1)
    mtx(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
end
a_diff = max(max((abs(a - mtx))))

as=sparse(a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating f_{ijk} and d_{ijk}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:(N+1)^2-1
    for j=1:(N+1)^2-1
        for k=1:(N+1)^2-1
            f_temp(j,k)=-imag1*trace(F{k+1}*(F{i+1}*F{j+1}-F{j+1}*F{i+1}));
            d_temp(j,k)=trace(F{k+1}*(F{i+1}*F{j+1}+F{j+1}*F{i+1}));
        end
    end
    fs{i}=sparse(f_temp);
    nzf(i)=nzmax(fs{i}); %% to test numel in sparse representation
    ds{i}=sparse(d_temp);
    nzd(i)=nzmax(ds{i}); %% to test numel
end
f_temp=0;
d_temp=0;

QEs=sparse(zeros((N+1)^2-1));
QUs=sparse(zeros((N+1)^2-1));
QJs=sparse(zeros((N+1)^2-1));
for i=1:(N+1)^2-1
    QEs=QEs+hE(i)*transpose(fs{i}); %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
    QJs=QJs+hJ(i)*transpose(fs{i}); %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
    QUs=QUs+hU(i)*transpose(fs{i}); %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
end

Qs = E0*QEs+U/N*QUs+J*QJs;
Qs_drv = QEs;

fn = sprintf('%s/Qs_base.txt', path_to_debug);
mtx_data = importdata(fn);
mtx = zeros((N+1)^2-1);
for st = 1:size(mtx_data, 1)
    mtx(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
end
Qs_diff = max(max((abs(Qs - mtx))))

fn = sprintf('%s/Qs_drv.txt', path_to_debug);
mtx_data = importdata(fn);
mtx = zeros((N+1)^2-1);
for st = 1:size(mtx_data, 1)
    mtx(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
end
Qs_drv_diff = max(max((abs(Qs_drv - mtx))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caclulating K_s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=zeros((N+1)^2-1,1);
for j=1:(N+1)^2-1
    k_temp=0;
    for i=1:(N+1)^2-1
        k_temp=k_temp+as(i,:)*fs{i}(:,j); %% to make a sum over middle index
    end
    K(j)=k_temp;
end
K=imag1/(N+1)*K;
Ks=sparse(K);

fn = sprintf('%s/Ks.txt', path_to_debug);
mtx_data = importdata(fn);
vec = zeros((N+1)^2-1, 1);
for st = 1:size(mtx_data, 1)
    vec(st) = mtx_data(st, 1) + imag1 * mtx_data(st, 2);
end
Ks_diff = max(max((abs(full(Ks) - vec))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caclulating R_{sm}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N.B. Matrix operations are used for acceleration. Result verified
% against direct multi-cycle calculation

R=zeros((N+1)^2-1);
for i=1:(N+1)^2-1
    for k=1:(N+1)^2-1
        R=R+as(i,k)*(transpose(fs{k})*(fs{i}+imag1*ds{i})...
            +transpose(fs{i})*conj(fs{k}+imag1*ds{k}));
    end
end
R=-R/4;
Rs=sparse(R);

Gs = E0*QEs+U/N*QUs+J*QJs + g/2*Rs;



A=A1-sqrt(-1)*A2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=0;
kk=kk+1;

H=(E0*HE+U/N*HU+J*HJ);

P=zeros((N+1)^2);
P=-sqrt(-1)*(kron(eye(N+1),H)-kron(transpose(H),eye(N+1)))+...
    g/2*(2*kron(eye(N+1),A)*kron(transpose(A'),eye(N+1))-...
    kron(transpose(A'*A),eye(N+1))-kron(eye(N+1),A'*A));
Ps=sparse(P);
P=0;

Pd=sparse(-sqrt(-1)*(kron(eye(N+1),HE)-kron(transpose(HE),eye(N+1))));



RhoFloquet=eye((N+1)^2 - 1);
diffs = zeros((N+1)^2 - 1, 1);
diffs_before = zeros((N+1)^2 - 1, 1);
diffs_cpp = zeros((N+1)^2 - 1, 1);
for jj=1:(N+1)^2 - 1
    
    t0=0;
    RhoF=RhoFloquet(:,jj);
    
    start_rho = F{1}/(N+1);
    for i=2:(N+1)^2
        start_rho=start_rho+RhoF(i-1)*F{i};
    end
    start_rho=full(start_rho);
    
    fn = sprintf('%s/rho_before_%d.txt', path_to_debug, jj - 1);
    mtx_data = importdata(fn);
    rho_mtx_before = zeros(N+1);
    for st = 1:size(mtx_data, 1)
        rho_mtx_before(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
    end
    diffs_before(jj) = max(max(abs(start_rho - rho_mtx_before)));
    if diffs_before(jj) > 1e-13
        diffs_before(jj) = max(max(abs(start_rho - transpose(rho_mtx_before))));
        if diffs_before(jj) < 1e-13
            start_rho = transpose(start_rho);
            RhoF = conj(RhoF);
        end
    end
        
        
    [t,Y]=ode45('sysQ',[t0:pi/w:t0+2*pi/w],RhoF);
    RhoFloquet(:,jj)=transpose(Y(length(t),:));
    rho_f = RhoFloquet(:,jj);
    Rho_fbasis_mtx=F{1}/(N+1);
    for i=2:(N+1)^2
        Rho_fbasis_mtx=Rho_fbasis_mtx+rho_f(i-1)*F{i};
    end
    Rho_fbasis_mtx=full(Rho_fbasis_mtx);
    
    start_rho_vec = zeros((N+1)^2, 1);
    for st_1 = 1:N+1
        for st_2 = 1:N+1
            index = (st_1-1) * (N+1) + st_2;
            start_rho_vec(index) = start_rho(st_1, st_2);
        end
    end

    [t,Y]=ode45('sysLind',[t0:pi/w:t0+2*pi/w],start_rho_vec);
    Rho_direct=transpose(Y(length(t),:));
    Rho_direct_mtx = zeros(N+1);
    for st_1 = 1:N+1
        for st_2 = 1:N+1
            Rho_direct_mtx(st_1, st_2) = Rho_direct((st_1-1) * (N+1) + st_2);
        end
    end
    
    diffs(jj) = max(max(abs(Rho_direct_mtx - Rho_fbasis_mtx)));
    
    fn = sprintf('%s/rho_after_%d.txt', path_to_debug, jj - 1);
    mtx_data = importdata(fn);
    rho_after_mtx_fb = zeros(N+1);
    for st = 1:size(mtx_data, 1)
        rho_after_mtx_fb(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
    end
    
    diffs_cpp(jj) = max(max(abs(Rho_direct_mtx - rho_after_mtx_fb)));

    t0=t(length(t));
end

toc

tic
eigv=eig(RhoFloquet);



