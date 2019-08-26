clear
global N J E0 A0 w U g HE HU HJ A phi0

N=50; % number of particles, system size = N+1
show=0; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);

tic

J=-1; % hopping constant
E0=-1.0; % bias
A0=-1.5; %driving amplitude
w=1; % driving frequency
T=2*pi/w;
phi0=pi*0; % initial phase
U=2.2/4; % on-site interaction
%Uarray=[5.:100.1:13]/1;
g=0.1/N; % gamma; 



%% NB: we will multiply the relevant coefficients by g/2 in the end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
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

Hp=((E0+A0)*HE+U/N*HU+J*HJ);

Hm=((E0-A0)*HE+U/N*HU+J*HJ);


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


%%%%%%%%%%%%%%%%%%%
% Creating matrices in the RHS
%%%%%%%%%%%%%%%%%%%
%-sqrt(-1)*(H*Rho-Rho*H)+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho)=0

Pp=zeros((N+1)^2);
Pp=-sqrt(-1)*(kron(eye(N+1),Hp)-kron(transpose(Hp),eye(N+1)))+...
g/2*(2*kron(eye(N+1),A)*kron(transpose(A'),eye(N+1))-...
kron(transpose(A'*A),eye(N+1))-kron(eye(N+1),A'*A));


Pm=zeros((N+1)^2);
Pm=-sqrt(-1)*(kron(eye(N+1),Hm)-kron(transpose(Hm),eye(N+1)))+...
g/2*(2*kron(eye(N+1),A)*kron(transpose(A'),eye(N+1))-...
kron(transpose(A'*A),eye(N+1))-kron(eye(N+1),A'*A));

RhoFloquet=expm(Pm*T/2)*expm(Pp*T/2);

[VFloquet,EFloquet]=eig(RhoFloquet);

%[BL,S,BR]=svd(P);
%ii=0;
%for i=1:(N+1)^2
%    if(abs(S(i,i))<1e-8)
%        ii=i;
%    end
%end

[fe index]=max((diag(real(EFloquet))))

Rho2=VFloquet(:,index);
for i=1:N+1
    Rho(:,i)=Rho2(1+(i-1)*(N+1):i*(N+1));
end

Rho=Rho/trace(Rho);

toc

dy=-sqrt(-1)*(Hp*Rho-Rho*Hp)+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho);
err=max(max(abs(dy)))

   Rho_plot=padarray(Rho,[1 1],'post');
   
   figure;
   [nx,ny]=meshgrid([0:N+1],[0:N+1]);
   h=pcolor(nx,ny,abs(Rho_plot));
   set(h,'EdgeColor','None')
   xlabel('m','FontSize',20)
   ylabel('n','FontSize',20)
   h1=gca;
   set(h1,'FontSize',16);
   title('N=20, J=1, E_0=0, \gamma=1, U=1         |\rho_{m,n}|','FontSize',16)
   colormap('hot')
   colorbar
   %print -dpng 'Figs_v1/abs_r_g1_U1.png'
   
   figure;
   %[nx,ny]=meshgrid([0:N],[0:N]);
   h=pcolor(nx,ny,angle(Rho_plot));
   set(h,'EdgeColor','None')
   xlabel('m','FontSize',20)
   ylabel('n','FontSize',20)
   h1=gca;
   set(h1,'FontSize',16);
   title('N=20, J=1, E_0=0, \gamma=1, U=1         arg(\rho_{m,n})','FontSize',16)
   colormap('hot')
   colorbar
   %print -dpng 'Figs_v1\angle_r_g1_U1.png'

   figure(6);
 hold on
 plot(diag(abs(Rho)),'o-b')
%plot(diff(diag(Rho)),'o-r')
%plot(diff(diag(Rho),2),'o-g')

  
 theta=linspace(0,pi,100);
 phi=linspace(-pi,pi,100);
 Hu=Husimi(theta,phi,Rho);  
 
 figure;
  [theta_grid,phi_grid]=meshgrid(theta,phi);
   h=pcolor(theta_grid,phi_grid,real(Hu'));
      set(h,'EdgeColor','None')
 
 figure;
 plot(diag(EFloquet),'o')
 
 
 break
 
 

[Vp,Ep,Wp]=eig(Pp);
[Vm,Em,Wm]=eig(Pm);

for k=1:(N+1)^2
    d=sum(Wp(:,k)'*Vp(:,k))^0.5;
    Vp(:,k)=Vp(:,k)/d;
    Wp(:,k)=Wp(:,k)/d;
    d=sum(Wm(:,k)'*Vm(:,k))^0.5;
    Vm(:,k)=Vm(:,k)/d;
    Wm(:,k)=Wm(:,k)/d;

end

RhoFloquet=zeros((N+1)^2);

for k=1:(N+1)^2
    Rho0=zeros((N+1)^2,1);
    Rho0(k,1)=1;

    C0=Wp'*Rho0;
        
    Rhop=zeros((N+1)^2,1);
    for kk=1:(N+1)^2
        Rhop=Rhop+C0(kk)*exp(Ep(kk,kk)*T/2)*Vp(:,kk);
    end

    Cp=Wm'*Rhop;
    
    for kk=1:(N+1)^2
        RhoFloquet(:,k)=RhoFloquet(:,k)+Cp(kk)*exp(Em(kk,kk)*T/2)*Vm(:,kk);
    end
    
end
