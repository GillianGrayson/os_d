clear;

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
Uarray=[0.07:0.005:0.15]/1;
g=0.1/N; % gamma;
eigvals=zeros(length(Uarray)*(N+1)^2,3);

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
    
end


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