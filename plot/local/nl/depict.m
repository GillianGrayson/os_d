clear
N=1000;

out=load('out.txt');

U=out(:,3);
A0=out(:,5);
t=out(:,7);
theta=out(:,8);
phi=out(:,9);

n1=N/2*(cos(theta)+1);

%figure;
%plot(U/4,n1,'.b','MarkerSize',0.25)
%plot(n1,'.b','MarkerSize',0.5)


%figure;
%plot(U,cos(theta)/2,'.','MarkerSize',1)

%figure;
%plot(U,sin(phi),'.','MarkerSize',1)


theta=mod(theta,2*pi);
phi=mod(phi,2*pi);

%figure;
%plot(theta,phi,'.','MarkerSize',0.5)

Ngrid=200;
LU=150;
Ln=length(U)/LU;
P=zeros(Ngrid,LU);
for k1=1:LU
    for k2=1:Ln
        id = floor( n1((k1-1)*Ln+k2) / N * (Ngrid-1)) + 1;
        P(id, k1)=P(id, k1)+1;
    end
    P(:,k1)=P(:,k1) / max(P(:,k1));
end

[theta_grid,phi_grid]=meshgrid(linspace(0,0.75,LU),linspace(0,N,Ngrid));
figure;
h=pcolor(theta_grid,phi_grid,P);
set(h,'EdgeColor','None')
colorbar


% Ngrid=200;
% P=zeros(Ngrid);
% for i=1:length(theta)
%     k_theta=floor(theta(i)/pi*(Ngrid-1))+1;
%     k_phi=floor(phi(i)/2/pi*(Ngrid-1))+1;
%     P(k_theta,k_phi)=P(k_theta,k_phi)+1;
% end
% P=P/sum(sum(P));
% 
% [theta_grid,phi_grid]=meshgrid(linspace(0,pi,Ngrid),linspace(0,2*pi,Ngrid));
% figure(5)
% h=pcolor(theta_grid,phi_grid,P');
% set(h,'EdgeColor','None')
% colorbar