clear;

N=3; % number of particles, system size = N+1
show=2; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);

tic

J=-1; % hopping constant
E0=-1; % on-site potential
U=2.2; % on-site interaction
g=0.1; % gamma;
%% NB: we will multiply the relevant coefficients by g/2 in the end

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
% Calculating H in F-basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=diag(E0*(2*[0:N]-N)+U/N*(2*[0:N]-N).^2); % diagonal
for i=1:N
    H(i+1,i)=J*sqrt((N-i+1)*(i));
    H(i,i+1)=J*sqrt((i)*(N-i+1));
end

H=H-trace(H)/(N+1)*eye(N+1);


Hs=sparse(H);
%% from here and further on we expand names of variables with "s"
%% to mark matrices in sparse representationc N^2 +2 N + 1

for i=1:(N+1)^2-1
    h(i)=trace(Hs*F{i+1});
end

if(show==1)
    H=H
    h=h
end

%break

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dissipator in F-basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1s=sparse(diag(2*[0:N]-N));
A2=zeros(N+1);
for i=1:N
    A2(i+1,i)=-sqrt((N-i+1)*(i));
    A2(i,i+1)=sqrt((i)*(N-i+1));
end
A2=imag1*A2;
A2s=sparse(A2);
A1 = full(A1s);
AA = A1 - imag1*A2;

EE = (eye(N+1)+ imag1* eye(N+1))* sqrt(1.0/ 2.0/(N+1))
tr = trace(EE);
EET = EE * ctranspose(AA)
R1 = AA * EET
R2 = EET * AA
R = R1 - R2
EET = EE * ctranspose(AA)
A2=0;

%A1s=sparse(diag([0:N]));
%A2s=sparse(diag([N:-1:0]));

if(show==1)
    full(A1s)
    full(A2s)
end

%break

for i=1:(N+1)^2-1
    a1(i)=trace(A1s*F{i+1});
    a2(i)=trace(A2s*F{i+1});
end

if(show==1)
    a1=a1
    a2=a2
end

for i=1:(N+1)^2-1
    for j=1:(N+1)^2-1
        a(i,j)=(a1(i)-imag1*a2(i))*conj(a1(j)-imag1*a2(j));
    end
end

%for i=1:(N+1)^2
%    for j=1:(N+1)^2
%        a(i,j)=a1(i)*a1(j)+a2(i)*a2(j);
%    end
%end


as=sparse(a);
a=0;

if(show==1)
    a=full(as)
    a=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating f_{ijk} and d_{ijk}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:(N+1)^2-1
    for j=1:(N+1)^2-1
        for k=1:(N+1)^2-1
            f_temp(j,k)=-imag1*trace(F{i+1}*(F{j+1}*F{k+1}-F{k+1}*F{j+1}));
            d_temp(j,k)=trace(F{i+1}*(F{j+1}*F{k+1}+F{k+1}*F{j+1}));
        end
    end
    fs{i}=sparse(f_temp);
    nzf(i)=nzmax(fs{i}); %% to test numel in sparse representation
    ds{i}=sparse(d_temp);
    nzd(i)=nzmax(ds{i}); %% to test numel
end
f_temp=0;
d_temp=0;

if(show==2)
    sumnzf=sum(nzf) %% actual numel of non-zero f_{ijk}
    sumnzd=sum(nzd) %% actual numel of non-zero d_{ijk}
    max_vol=(N+1)^6 %% numel of f_{ijk} and d_{ijk} including zeros
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caclulating Q_{sm}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qs=sparse(zeros((N+1)^2-1));
for i=1:(N+1)^2-1
    Qs=Qs+h(i)*fs{i}; %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
end

%if(show==1)
Q_full=full(Qs)
Q_full=0;
%end

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
K=g/N*imag1/(N+1)*K;
Ks=sparse(K);
K=0;

%   if(show==1)
K_full=full(Ks)
K_full=0;
%   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caclulating R_{sm}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N.B. Matrix operations are used for acceleration. Result verified
% against direct multi-cycle calculation


R=zeros((N+1)^2-1);
for i=1:(N+1)^2-1
    for k=1:(N+1)^2-1
        R=R+as(i,k)*(fs{k}.'*(fs{i}+imag1*ds{i})+fs{i}.'*conj(fs{k}+imag1*ds{k}));
    end
end
R=-R/4*g/N;
Rs=sparse(R);
R=0;

%   if(show==1)
R=full(Rs)
R=0;
%  end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving for a stationary point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gs=Qs+Rs;
if(max(max(abs(imag(Gs))))<1e-15)
    Gs=real(Gs)
end

if(max(abs(imag(Ks)))<1e-15)
    Ks=real(Ks)
end


%G_to_solve=full(Gs(2:(N+1)^2,2:(N+1)^2));
%K_to_solve=full(Ks(2:(N+1)^2));
%det(G_to_solve)

%if(detG~=0)
%RhoF=linsolve(G_to_solve,-K_to_solve)
RhoF=linsolve(full(Gs),-full(Ks))

%else
%% Two of many solutions for a consistent det=0 systems
%if(rank([full(Gs),-full(Ks)])==rank(full(Gs)))
%    RhoF1=-Ks\Gs
%    RhoF2=-pinv(full(Gs))*Ks
%
%end
%end

Rho=F{1}/(N+1);
for i=2:(N+1)^2
    Rho=Rho+RhoF(i-1)*F{i};
end

Rho=full(Rho)

trace(Rho)
trace(Rho^2)

toc

figure;
[nx,ny]=meshgrid([0:N],[0:N]);
h=pcolor(nx,ny,abs(Rho));
set(h,'EdgeColor','None')
xlabel('m')
ylabel('n')
title('abs(\rho_{m,n})')
colorbar

figure;
[nx,ny]=meshgrid([0:N],[0:N]);
h=pcolor(nx,ny,angle(Rho));
set(h,'EdgeColor','None')
xlabel('m')
ylabel('n')
title('arg(\rho_{m,n})')
colorbar