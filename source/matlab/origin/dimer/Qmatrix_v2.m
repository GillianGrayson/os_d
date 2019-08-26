clear;

N=2; % number of particles, system size = N+1
show=0; % if = 1, generated matrices, coefficients, etc. are displayed
imag1=sqrt(-1);

tic

%J=1; % hopping constant
%E0=0; % on-site potential
%U=1; % on-site interaction
%g=10; % gamma; 
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
HE=zeros(N+1);
HJ=zeros(N+1);
HU=zeros(N+1);

HE=diag(2*[0:N]-N); % diagonal
HU=diag(1/N*(2*[0:N]-N).^2); % diagonal
for i=1:N
    HJ(i+1,i)=sqrt((N-i+1)*(i));
    HJ(i,i+1)=sqrt((i)*(N-i+1));
end

HE=HE-trace(HE)/(N+1)*eye(N+1);
HU=HU-trace(HU)/(N+1)*eye(N+1);
HJ=HJ-trace(HJ)/(N+1)*eye(N+1);


HEs=sparse(HE); 
HUs=sparse(HU); 
HJs=sparse(HJ); 


%% from here and further on we expand names of variables with "s" 
%% to mark matrices in sparse representation
  
for i=1:(N+1)^2-1
    hE(i)=trace(HEs*F{i+1});
    hU(i)=trace(HUs*F{i+1});
    hJ(i)=trace(HJs*F{i+1});
end

if(show==1)
    HE=HE
    hE=hE
    HU=HU
    hU=hU
    HJ=HJ
    hJ=hJ
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
 %% !!!!!!!!
               f_temp(j,k)=-imag1*trace(F{k+1}*(F{i+1}*F{j+1}-F{j+1}*F{i+1}));
               d_temp(j,k)=trace(F{k+1}*(F{i+1}*F{j+1}+F{j+1}*F{i+1}));
 %% ????????
     % f_temp(j,k)=-imag1*trace(F{i+1}*(F{j+1}*F{k+1}-F{k+1}*F{j+1}));
      %         d_temp(j,k)=trace(F{i+1}*(F{j+1}*F{k+1}+F{k+1}*F{j+1}));
          

           end
       end
       fs{i}=sparse(f_temp);
       nzf(i)=nzmax(fs{i}); %% to test numel in sparse representation
       ds{i}=sparse(d_temp); 
       nzd(i)=nzmax(ds{i}); %% to test numel
   end
   f_temp=0;
   d_temp=0;
   
   if(show==1)
   sumnzf=sum(nzf) %% actual numel of non-zero f_{ijk}
   sumnzd=sum(nzd) %% actual numel of non-zero d_{ijk}
   max_vol=(N+1)^6 %% numel of f_{ijk} and d_{ijk} including zeros
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Caclulating Q_{sm}
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   QEs=sparse(zeros((N+1)^2-1));
   QUs=sparse(zeros((N+1)^2-1));
   QJs=sparse(zeros((N+1)^2-1));
   for i=1:(N+1)^2-1
 %%%% !!!!!!!!!!!
      QEs=QEs+hE(i)*transpose(fs{i}); %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
      QJs=QJs+hJ(i)*transpose(fs{i}); %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
      QUs=QUs+hU(i)*transpose(fs{i}); %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
 %%%%% ??????????
      %QEs=QEs+hE(i)*fs{i}; %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
      %QJs=QJs+hJ(i)*fs{i}; %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2
      %QUs=QUs+hU(i)*fs{i}; %% fs{i} is a sparse matrix of f_{i,j,k}, i fixed, j,k=1...(N+1)^2

   
   end
   
   if(show==1)
       QE_full=full(QEs)
       QE_full=0;
       QU_full=full(QUs)
       QU_full=0;
       QJ_full=full(QJs)
       QJ_full=0;
   end
   
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
   %%%%%% !!!!!!!!!
   K=imag1/(N+1)*K;
   %%%%%% ?????????
   Ks=sparse(K);
   K=0;
   
   if(show==1)
       K_full=full(Ks)
       K_full=0;
   end
   
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
   %%%%%%%%% !!!!!!!!!!
   R=-R/4;
   %%%%%%%%% ?????????
   Rs=sparse(R);
   R=0;
   
   if(show==1)
       R=full(Rs)
       R=0;
   end
  
   toc
   
       QE_full=real(full(QEs));
       save QE.txt QE_full -ascii -double
       QE_full=0;
  
       QU_full=real(full(QUs));
       save QU.txt QU_full -ascii -double
       QU_full=0;
  
       QJ_full=real(full(QJs));
       save QJ.txt QJ_full -ascii -double
       QJ_full=0;
  
        K_full=real(full(Ks));
        save K.txt K_full -ascii -double
        K_full=0;
 
        R=real(full(Rs));
        save R.txt R -ascii -double
        R=0;