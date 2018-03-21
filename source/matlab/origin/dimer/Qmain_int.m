clear;

global N J E0 A0 w U g QEs QUs QJs Rs Ks phi

N=20; % number of particles, system size = N+1
show=0; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);

tic

J=1; % hopping constant
E0=0; % bias
A0=10; %driving amplitude
w=5; % driving frequency
phi=pi*0; % initial phase
U=1.; % on-site interaction
g=1; % gamma;

N_it_pre=100; % number of transient periods
N_it=10; % number of plotted time points

%% NB: we will multiply the relevant coefficients by g/2 in the end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QE=load('QE.txt');
QEs=sparse(QE);
QE=0;
QU=load('QU.txt');
QUs=sparse(QU);
QU=0;
QJ=load('QJ.txt');
QJs=sparse(QJ);
QJ=0;

R=load('R.txt');
Rs=sparse(R);
R=0;

K=load('K.txt');
Ks=sparse(K);
K=0;


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

if(show==1)
    for i=1:(N+1)^2
        m_show=full(F{i})
    end
    
    for i=1:(N+1)^2
        trace_kk(i)=trace(F{i});
    end
    
    for i=1:(N+1)^2
        for j=1:(N+1)^2
            trace_kj(i,j)=trace(F{i}*F{j});
        end
    end
    
    trace_kk=trace_kk
    trace_kj=trace_kj
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gs=(E0*QEs+U*QUs+J*QJs)+g*Rs;
%Ks=g*Ks;

%if(max(max(abs(imag(Gs))))<1e-15)
%    Gs=real(Gs);
%end

%if(max(abs(imag(Ks)))<1e-15)
%    Ks=real(Ks);
%end


%       RhoF=linsolve(full(Gs),-full(Ks));

%    eVal=eig(full(Gs));
RhoF=rand((N+1)^2-1,1);
t0=0;
for jj=1:N_it_pre
    
    [t,Y]=ode45('sysQ',[t0:pi/w:t0+2*pi/w],RhoF);
    RhoF=Y(length(t),:)';
    t0=t(length(t));
end

for jj=N_it_pre:N_it_pre+N_it
    
    [t,Y]=ode45('sysQ',[t0:pi/w/N_it:t0+pi/w*2/N_it],RhoF);
    RhoF=Y(length(t),:)';
    t0=t(length(t));
    
    
    Rho=F{1}/(N+1);
    for i=2:(N+1)^2
        Rho=Rho+RhoF(i-1)*F{i};
    end
    
    Rho=full(Rho);
    
    trace(Rho);
    trace(Rho^2);
    
    
    Rho_plot=padarray(Rho,[1 1],'post');
    
    figure;
    
    subplot(2,2,1)
    
    [nx,ny]=meshgrid([0:N+1],[0:N+1]);
    h=pcolor(nx,ny,(abs(Rho_plot)));
    set(h,'EdgeColor','None')
    xlabel('m')
    ylabel('n')
    title('abs(\rho_{m,n})')
    colormap('hot')
    colorbar
    
    subplot(2,2,2)
    
    [nx,ny]=meshgrid([0:N+1],[0:N+1]);
    h=pcolor(nx,ny,(angle(Rho_plot)));
    set(h,'EdgeColor','None')
    xlabel('m')
    ylabel('n')
    title('abs(\rho_{m,n})')
    colormap('hot')
    colorbar
    
    % figure;
    % %[nx,ny]=meshgrid([0:N],[0:N]);
    % h=pcolor(nx,ny,angle(Rho_plot));
    % set(h,'EdgeColor','None')
    % xlabel('m')
    % ylabel('n')
    % title('arg(\rho_{m,n})')
    % colormap('hot')
    % colorbar
    
end

toc

