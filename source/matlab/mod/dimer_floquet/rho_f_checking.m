clear;

global N J E0 A0 w U g HE HU HJ A phi0 Ps Pd
imag1=sqrt(-1);


N=10; % number of particles, system size = N+1

path_to_debug = 'E:/Work/os_d/source/cpp/CQdiss_os_dimer/CQdiss_fbasis';

tic

J=-1; % hopping constant
E0=0.; % bias
A0=-3.4; %driving amplitude
w=1; % driving frequency
phi0=pi*0; % initial phase
U=0.5; % on-site interaction
g=0.1/N; % gamma;
eigvals=zeros((N+1)^2,3);

N_it_pre=100; % number of transient periods
N_it_post=100; % number of plotted periods
N_per=10;


%%%%%%%%%%%%%%%%%%%%%%%% F-basis matrixes %%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating H
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dissipator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1=diag(2*[0:N]-N);
A2=zeros(N+1);
for i=1:N
    A2(i+1,i)=-sqrt((N-i+1)*(i));
    A2(i,i+1)=sqrt((i)*(N-i+1));
end
A2=imag1*A2;

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

toc
tic

rho_f_sim = (N+1)^2 - 1;
RhoFloquet = zeros(rho_f_sim);
diffs = zeros(rho_f_sim, 1);

for jj=1:rho_f_sim
    
    t0 = 0;
    
    fn = sprintf('%s/rho_before_%d.txt', path_to_debug, jj - 1);
    mtx_data = importdata(fn);
    rho_mtx_before = zeros(N+1);
    for st = 1:size(mtx_data, 1)
        rho_mtx_before(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
    end
    rho_before = zeros((N+1)^2, 1);
    for st_1 = 1:N+1
        for st_2 = 1:N+1
            index = (st_1-1) * (N+1) + st_2;
            rho_before(index) = rho_mtx_before(st_1, st_2);
        end
    end
    
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t,Y]=ode45('sysLind',[t0:pi/w:t0+2*pi/w],rho_before, options);
    rho_after = transpose(Y(length(t),:));
    
    %RhoFloquet(:, jj) = rho_after;
    
    rho_after_mtx = zeros(N+1);
    for st_1 = 1:N+1
        for st_2 = 1:N+1
            index = (st_1-1) * (N+1) + st_2;
            rho_after_mtx(st_1, st_2) = rho_after(index);
        end
    end
    
    fn = sprintf('%s/rho_after_%d.txt', path_to_debug, jj - 1);
    mtx_data = importdata(fn);
    rho_after_mtx_fb = zeros(N+1);
    for st = 1:size(mtx_data, 1)
        rho_after_mtx_fb(mtx_data(st, 1), mtx_data(st, 2)) = mtx_data(st, 3) + imag1 * mtx_data(st, 4);
    end
    
    abs_diff = max(max(abs(rho_after_mtx - rho_after_mtx_fb)));
    diffs(jj) =  abs_diff;
        
    t0 = t(length(t));
end

max(diffs)

toc

tic
eigv=eig(RhoFloquet);
toc
tic
[RhoF,S]=eigs(RhoFloquet,1,'lm');
toc
for i=1:N+1
    Rho(:,i)=RhoF(1+(i-1)*(N+1):i*(N+1),1);
end
Rho=Rho/trace(Rho);
Purity=trace(Rho^2)
Negativity=0.5*(sum(sum(abs(Rho)))-trace(Rho))


figure;
plot(eigv,'og')
hold on
plot(cos([0:0.01:2*pi]),sin([0:0.01:2*pi]),'-')

temp=sort(abs(eigv));


eigmax(kk)=temp(length(temp)-1);
eigmin(kk)=temp(1);

eigvals(1+(kk-1)*(N+1)^2:kk*(N+1)^2,1)=U*ones((N+1)^2,1);
eigvals(1+(kk-1)*(N+1)^2:kk*(N+1)^2,2)=real(eigv);
eigvals(1+(kk-1)*(N+1)^2:kk*(N+1)^2,3)=imag(eigv);
delta=1-eigmax(kk)



figure;

plot(Uarray,eigmin,'s-b')
hold on
plot(Uarray,eigmax,'s-r')

%Rho2=Rho2/N_it_post;
%RhoD=RhoD/N_it_post;
%eVmax=eVmax/N_it_post;
tic

theta=linspace(0,pi,50);
phi=linspace(0,2*pi,50);
Hu=Husimi(theta,phi,Rho);

toc

figure;
[theta_grid,phi_grid]=meshgrid(theta,phi);
h=pcolor(theta_grid,phi_grid,real(Hu'));
set(h,'EdgeColor','None')




toc


save eigvals.txt eigvals -ascii