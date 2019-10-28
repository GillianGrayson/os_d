clear all;

global N J E0 A0 w U g HE HU HJ A phi0 Ps Pd

N=50; % number of particles, system size = N+1
show=0; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);

tic

J=-1; % hopping constant
E0=0.; % bias
A0=-3.4; %driving amplitude
w=1; % driving frequency
phi0=pi*0; % initial phase
%U=0.7; % on-site interaction
Uarray=[0.12];
g=0.1/N; % gamma;
eigvals=zeros((N+1)^2, length(Uarray));

N_it_pre=100; % number of transient periods
N_it_post=100; % number of plotted periods
N_per=10;

%RhoD=zeros(N+1,length(N_it_post));
%eVmax=zeros(length(N_it_post));
%Rho2=zeros(length(N_it_post));


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
for U=Uarray
    kk=kk+1
    
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
    
    RhoFloquet=eye((N+1)^2);
    
    
    
    for jj=1:(N+1)^2
        
        t0=0;
        RhoF=RhoFloquet(:,jj);
        
        [t,Y]=ode45('sysLind',[t0:pi/w:t0+2*pi/w],RhoF);
        RhoFloquet(:,jj)=transpose(Y(length(t),:));
        t0=t(length(t));
    end
    
    toc
    
    tic
    eigv=eig(RhoFloquet);
    eigvals(:, kk) = eigv;
    toc
    
end
