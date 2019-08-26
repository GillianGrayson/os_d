clear;

global N J E0 A0 w U g HE HU HJ A phi0

N=20; % number of particles, system size = N+1
show=0; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);

tic

J=-1; % hopping constant
E0=-1.; % bias
A0=-1.5; %driving amplitude
w=1; % driving frequency
phi0=pi*0; % initial phase
%U=0.7; % on-site interaction
Uarray=[0.5:100.1:30]/4;
g=0.1/N; % gamma; 

N_it_pre=100; % number of transient periods
N_it_post=100; % number of plotted periods
N_per=10; 

RhoD=zeros(N+1,length(N_it_post));
eVmax=zeros(length(N_it_post));
Rho2=zeros(length(N_it_post));


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

H=(E0*HE+U/N*HU+J*HJ);

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
       
   %psi0=rand(N+1,1)+sqrt(-1)*rand(N+1,1);
   psi0=ones(N+1,1)+sqrt(-1)*ones(N+1,1);
   psi0=psi0/sqrt(sum(abs(psi0.^2)));
   Rho0=psi0*psi0';
   
   for i=1:N+1
    RhoF(1+(i-1)*(N+1):i*(N+1),1)=Rho0(:,i);
   end
   
   t0=0;
    for jj=1:N_it_pre
       
       [t,Y]=ode45('sysLind',[t0:pi/w:t0+2*pi/w],RhoF);
       RhoF=transpose(Y(length(t),:));
       t0=t(length(t));
    end
   
   for jj=1:N_per*N_it_post
       
       [t,Y]=ode45('sysLind',[t0:pi/w/N_it_post:t0+2*pi/w/N_it_post],RhoF);
       RhoF=transpose(Y(length(t),:));
       t0=t(length(t));
        
  
  
for i=1:N+1
    Rho(:,i)=RhoF(1+(i-1)*(N+1):i*(N+1),1);
end

   trace(Rho);
   trace(Rho^2);
   
   Rho2(jj)=trace(Rho^2);
   
   RhoD(:,jj)=diag(Rho);
   %eVmax(1,kk)=eVmax(1,kk)+max(real(eVal));
   
   
   end
  
   end
   
   %Rho2=Rho2/N_it_post;
   %RhoD=RhoD/N_it_post;
   %eVmax=eVmax/N_it_post;
   toc
       
   
   Rho_plot=padarray(RhoD,[1 1],'post');
  
   Narray=[1:N_it_post*N_per];
   
   figure;
   [nx,ny]=meshgrid([Narray 2*Narray(length(Narray))-Narray(length(Narray)-1)],[0:N+1]);
   h=pcolor(nx,ny,(abs(Rho_plot)));
   set(h,'EdgeColor','None')
 %mesh(nx,ny,abs(Rho_plot));
 axis tight 
 
  xlabel('U','FontSize',20)
   ylabel('n','FontSize',20)
   title('N=20, J=1.0, E_0=0.05, \gamma=1,       |\rho_{n,n}|','FontSize',16)
   h1=gca;
   set(h1,'FontSize',16);
   colormap('hot')
   colorbar
  
   %print -dpng 'Figs_v1\rnn_U_J1_g1_E005.png'
   
   figure;
   title('N=20, J=1.0, E_0=0.05, \gamma=1','FontSize',16)
   subplot(2,1,1)
  % plot(U_array,log10(-eVmax(1,:)),'o-')
   plot(Uarray,eVmax(1,:),'o-')
  
  box on
   xlabel('U','FontSize',20)
   ylabel('log_{10}[-max Re(eig)]','FontSize',16)
   h1=gca;
   set(h1,'FontSize',16);
   
   subplot(2,1,2)
   plot(Uarray,Rho2,'o-')
   
   box on
   xlabel('U','FontSize',20)
   ylabel('Tr(\rho^2)','FontSize',20)
   h1=gca;
   set(h1,'FontSize',16);
   ylim([0 1]);
   %print -dpng 'Figs_v1\Lambda_U_J1_g1_E005.png'
   
   theta=linspace(0,pi,100);
 phi=linspace(-pi,pi,100);
 Hu=Husimi(theta,phi,Rho);  
 
 figure;
   [theta_grid,phi_grid]=meshgrid(theta,phi);
   h=pcolor(theta_grid,phi_grid,real(Hu'));
      set(h,'EdgeColor','None')

      
      
      
      
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
   
   
   
  break 
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Check
   %%%%%%%%%%%%%%%%%%%%%%%%%
   
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

ft=sin(w*t0+phi0);
H=((E0+1*A0*ft)*HE+U/N*HU+J*HJ);

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


dy=-sqrt(-1)*(H*Rho-Rho*H)+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho);
err=max(max(abs(dy)))

%%%%%%%%%%
dRhoF=sysLind(t0,RhoF);
   dRho=zeros(N+1);
   for i=2:(N+1)^2
       dRho=dRho+dRhoF(i-1)*F{i};
   end
   
derror=max(max(abs(dy-dRho)))