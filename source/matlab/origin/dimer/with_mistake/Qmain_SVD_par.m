clear
global N J E0 A0 w U g HE HU HJ A phi0

N=10; % number of particles, system size = N+1
show=0; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);


J=0.02; % hopping constant
E0=-0.0; % bias
A0=-0; %driving amplitude
w=1; % driving frequency
phi0=pi*0; % initial phase
U=0.7; % on-site interaction
U_array=0.:0.01:0.4;
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

k=0;
for U=U_array
    k=k+1;
    %%%%%%%%%%%%%%%%%%%
    % Creating matrix P*Rho=0
    %%%%%%%%%%%%%%%%%%%
    %-sqrt(-1)*(H*Rho-Rho*H)+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho)=0
    
    
    H=(E0*HE+U/N*HU+J*HJ);
    
    
    P=zeros((N+1)^2);
    P=-sqrt(-1)*(kron(eye(N+1),H)-kron(transpose(H),eye(N+1)))+...
        g/2*(2*kron(eye(N+1),A)*kron(transpose(A'),eye(N+1))-...
        kron(transpose(A'*A),eye(N+1))-kron(eye(N+1),A'*A));
    
    %[BL,S,BR]=svd(P);
    %ii=0;
    %for i=1:(N+1)^2
    %    if(abs(S(i,i))<1e-8)
    %        ii=i;
    %    end
    %end
    
    S0=sort(real(eig(P)),'descend');
    S1(k)=real(S0(1));
    S2(k)=real(S0(2));
    S3(k)=real(S0(3));
    
    Ps=sparse(P);
    P=0;
    [BR,S]=eigs(Ps,1,'sm');
    S=0;
    Ps=0;
    BR1=0;
    
    Rho2=BR(:,1);%(:,ii);
    BR=0;
    
    for i=1:N+1
        Rho(:,i)=Rho2(1+(i-1)*(N+1):i*(N+1));
    end
    Rho2=0;
    
    Rho=Rho/trace(Rho);
    
    RhoP(k)=trace(Rho^2);
    RhoD(:,k)=diag(Rho);
    
    toc
    
    dy=-sqrt(-1)*(H*Rho-Rho*H)+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho);
    err=max(max(abs(dy)))
    
end


Rho_plot=padarray(RhoD,[1 1],'post');

figure;
[nx,ny]=meshgrid([U_array 2*U_array(length(U_array))-U_array(length(U_array)-1)],[0:N+1]);
%h=pcolor(nx,ny,log10(abs(Rho_plot)+1e-3));
h=pcolor(nx,ny,(abs(Rho_plot)));
set(h,'EdgeColor','None')
%mesh(nx,ny,abs(Rho_plot));
axis tight

xlabel('U','FontSize',20)
ylabel('n','FontSize',20)
title('N=50, J=1.0, E_0=0.05, \gamma=1','FontSize',16)
h1=gca;
set(h1,'FontSize',16);
colormap('hot')
colorbar


figure;
[nx,ny]=meshgrid([U_array 2*U_array(length(U_array))-U_array(length(U_array)-1)],[0:N+1]);
h=pcolor(nx,ny,log10(abs(Rho_plot)+1e-3));
%h=pcolor(nx,ny,(abs(Rho_plot)));
set(h,'EdgeColor','None')
%mesh(nx,ny,abs(Rho_plot));
axis tight

xlabel('U','FontSize',20)
ylabel('n','FontSize',20)
title('N=50, J=1.0, E_0=0.05, \gamma=1','FontSize',16)
h1=gca;
set(h1,'FontSize',16);
colormap('hot')
colorbar



figure;
title('N=50, J=1.0, E_0=0.05, \gamma=1','FontSize',16)
subplot(2,1,1)
%plot(U_array,log10(-eVmax(1,:)),'o-')
plot(U_array,S1,'o-')
hold on
plot(U_array,S2,'o-r')
plot(U_array,S3,'o-g')


box on
xlabel('U','FontSize',20)
ylabel('\lambda_2','FontSize',16)
h1=gca;
set(h1,'FontSize',16);

subplot(2,1,2)
plot(U_array,RhoP,'o-')

box on
xlabel('U','FontSize',20)
ylabel('Tr(\rho^2)','FontSize',20)
h1=gca;
set(h1,'FontSize',16);
ylim([0 1]);

figure(6);
hold on
plot(diag(abs(Rho)),'o-b')
%plot(diff(diag(Rho)),'o-r')
%plot(diff(diag(Rho),2),'o-g')


%theta=linspace(0,pi,100);
%phi=linspace(0,2*pi,100);
%Hu=Husimi(theta,phi,Rho);

%figure;
%  [theta_grid,phi_grid]=meshgrid(theta,phi);
%  h=pcolor(theta_grid,phi_grid,real(Hu'));
%     set(h,'EdgeColor','None')