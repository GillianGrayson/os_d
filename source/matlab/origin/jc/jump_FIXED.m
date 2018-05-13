clear all; clc;

alfa=5;

xi=1;
gamma=0.25*xi;
gamma1=0.25/alfa;
gamma2=0.25/alfa;
n_sr=0;
%L=200; %?????????? ??????

%?????? ???????? ???????? (?_) ? ??????????? (?)
N=20; %?????? ?????? a ? a_
a=zeros(N,N);
for n=1:(N-1)
    a(n,n+1)=sqrt(n);
end
a_=a';
%????????????? ?????????
A1=sqrt(n_sr+1)*a;
A2=sqrt(n_sr)*a_;
%??????????? ???????????? H(t)=H1+F(t)*H2
H1=1/2*(1/(alfa*alfa*alfa))*a_*a_*a*a-(i/2)*(gamma1*A1'*A1+gamma2*A2'*A2);
H2=i*(a_-a);

%????????? ??? F(t)
f0=3.2;
F0=f0*xi;
%tau1=0.98;
%T=1.98;
tau1=0.98*alfa;
T=1.98*alfa;
t1=tau1/xi; %????? ??????????, ????? F ?? ????? 0
t2=T/xi-t1; %????? ??????????, ????? F ????? 0



K=200; %?????????? ????????, ??????? ????????????? (transient)

M_int=200; %???-?? ?????????? ??? ???????????
Prob=zeros(M_int,M_int);


t=0; %????????? ?????
c=zeros(N,1); c(1,1)=1; %????????? ?????? ?????????
t_crit=t1;

%********************************************************************** ???? ?????????? (? ???????? ?????????????)
NNNt=500; % ????? ????? ?? ?????? ? ??? ?????????????? ??????????
TT=T;
h=TT/NNNt; %??? ?????????????? ??? ??????????

NNNt1=int32(NNNt*t1/T);
NNNt2=NNNt-NNNt1;

G2=expm(-i*(H1)*h);
G1=expm(-i*(H1+F0*H2)*h);

t=0;
eta_sqrt=sqrt(rand);

for klkl = 1:K
    %****************************************************************************************************+
    
    %***********************************************************************************************
    for  kst = 1:NNNt1
        
        c=G1*c;
        norm_c=norm(c);
        
        
        if (norm_c < eta_sqrt)
            
            cs1=A1*c;
            c=cs1/norm(cs1);
            norm_c=1;
            eta_sqrt=sqrt(rand);
            
            
        end
        
        
    end
    
    
    
    for  kst = 1:NNNt2
        
        c=G2*c;
        norm_c=norm(c);
        
        
        if (norm_c < eta_sqrt)
            
            cs1=A1*c;
            c=cs1/norm(cs1);
            norm_c=1;
            eta_sqrt=sqrt(rand);
            
            
        end
        
        
    end
    
    
    
    %*****************************************************************************************************
    
end
%?????? ?????????????
%*****************************************************************************************************
NNN=400; % ????? ????? ?? ?????? ??? ????????????? ??????????
h=TT/NNN;
t=0;

NNN1=int32(NNN*t1/T);
NNN2=NNN-NNN1;

%GF = cell(NHier, 1)
%GL = cell(NHier, 1)

hs=h;
for  kss = 1:4
    GL{kss}=expm(-i*(H1)*hs);
    GF{kss}=expm(-i*(H1+F0*H2)*hs);
    hs=hs/10;
end


XXMax=12;
YYMax=12;



countR=1

eta_sqrt=sqrt(rand);
norm_c=norm(c);

for klls=1:50% ?????? ??????? ????... ????? ???????, ???? ????, ? ???????, ??????? ?????? ?????
    
    for klkl=1:1000 % ?????  ?????????? (?????????????? ?? ??????), ????? ??????? ??????????? ???? ? ??????? ???????????
        %****************************************************************************************************+
        click=0;
        flag=0;
        t=0;
        
        %************************************************************************** start of the first part of the period (force is acting)
        while flag < 1
            
            tu=0;
            
            c1=GF{1}*c;
            norm_c=norm(c1);
            
            
            if (norm_c < eta_sqrt)
                
                click=1;
                
                tu=0;
                
                %**************************************************************** start of the 2th level
                c2=c;
                
                while click > 0
                    
                    
                    for fer1=1:10                      %6
                        cf2=GF{2}*c2;
                        norm_c2=norm(cf2);
                        
                        tu=tu+h/10;
                        
                        
                        
                        if (norm_c2 < eta_sqrt)     %5
                            
                            
                            %************************************************************** start of the 3th level
                            c3=c2;
                            
                            while click > 0
                                
                                
                                for fer3=1:10                    %4
                                    cf3=GF{3}*c3;
                                    norm_c3=norm(cf3);
                                    
                                    
                                    tu=tu+h/100;
                                    
                                    if (norm_c3 < eta_sqrt)     %3
                                        
                                        
                                        %******************************************************* core of the hiearchy
                                        c4=c3;
                                        
                                        while click > 0
                                            
                                            for fer4=1:10                 %2
                                                
                                                cf4=GF{4}*c4;
                                                norm_c4=norm(cf4);
                                                
                                                tu=tu+h/1000;
                                                
                                                if (norm_c4 < eta_sqrt)  %1
                                                    
                                                    cs1=A1*c4;
                                                    cf4=cs1/norm(cs1);
                                                    eta_sqrt=sqrt(rand);
                                                    norm_c=1;
                                                    norm_c1=1;
                                                    norm_c2=1;
                                                    norm_c3=1;
                                                    norm_c4=1;
                                                    
                                                    click=0;
                                                end                         %1
                                                
                                                c4=cf4;
                                            end
                                        end                   %2
                                        %************************************************************** end of the core
                                        cf3=c4;
                                        
                                        
                                    end                        %3
                                    
                                    
                                    c3=cf3;
                                end
                            end        %4
                            %************************************************************************ end of the 3-th level
                            cf2=c3;
                        end
                        
                        
                        c2=cf2;
                        
                    end
                end
                %*************************************************************************** end of the 2-th level
                
                c1=c4;
                
            end
            
            
            if tu == 0
                tu=h;
            end
            
            tcur=t+tu;
            
            if tcur < t1
                flag=0;
                c=c1;
                t=t+tu;
            else
                dt=t1-t;
                t=t+dt;
                PFinal=expm(-i*(H1+F0*H2)*dt);
                flag=1;
                c=PFinal*c1;
            end
            
            
        end
        %************************************************************************************* end of first part of the period (force is acting)
        
        
        
        
        
        %************************************************************************************** beggining of the second part
        
        
        click=0;
        flag=0;
        t=0;
        
        %************************************************************************** start of the second part of the period (force is NOT acting)
        while flag < 1
            
            tu=0;
            
            c1=GL{1}*c;
            norm_c=norm(c1);
            
            
            if (norm_c < eta_sqrt)
                
                click=1;
                
                tu=0;
                
                %**************************************************************** start of the 2th level
                c2=c;
                
                while click > 0
                    
                    
                    for fer1=1:10                      %6
                        cf2=GL{2}*c2;
                        norm_c2=norm(cf2);
                        
                        tu=tu+h/10;
                        
                        
                        
                        if (norm_c2 < eta_sqrt)     %5
                            
                            
                            %************************************************************** start of the 3th level
                            c3=c2;
                            
                            while click > 0
                                
                                
                                for fer3=1:10                    %4
                                    cf3=GL{3}*c3;
                                    norm_c3=norm(cf3);
                                    
                                    
                                    tu=tu+h/100;
                                    
                                    if (norm_c3 < eta_sqrt)     %3
                                        
                                        
                                        %******************************************************* core of the hiearchy
                                        c4=c3;
                                        
                                        while click > 0
                                            
                                            for fer4=1:10                 %2
                                                
                                                cf4=GL{4}*c4;
                                                norm_c4=norm(cf4);
                                                
                                                tu=tu+h/1000;
                                                
                                                if (norm_c4 < eta_sqrt)  %1
                                                    
                                                    cs1=A1*c4;
                                                    cf4=cs1/norm(cs1);
                                                    eta_sqrt=sqrt(rand);
                                                    norm_c=1;
                                                    norm_c1=1;
                                                    norm_c2=1;
                                                    norm_c3=1;
                                                    norm_c4=1;
                                                    
                                                    click=0;
                                                end                         %1
                                                
                                                c4=cf4;
                                            end
                                        end                   %2
                                        %************************************************************** end of the core
                                        cf3=c4;
                                        
                                        
                                    end                        %3
                                    
                                    
                                    c3=cf3;
                                end
                            end        %4
                            %************************************************************************ end of the 3-th level
                            cf2=c3;
                        end
                        
                        
                        c2=cf2;
                        
                    end
                end
                %*************************************************************************** end of the 2-th level
                
                c1=c4;
                
            end
            
            
            if tu == 0
                tu=h;
            end
            
            tcur=t+tu;
            
            if tcur < t2
                flag=0;
                c=c1;
                t=t+tu;
            else
                dt=t2-t;
                PFinal=expm(-i*H1*dt);
                flag=1;
                c=PFinal*c1;
            end
            
            
        end
        
        %************************************************************************************* end of second part of the period (force is NOT acting)
        
        %************************************************************************************************
        
        zz=c/norm(c);
        xx = real(zz'*a*zz);
        yy = imag(zz'*a*zz);
        %  ?????? ??????? ?????? ??? ???????????
        if (abs(xx)<(XXMax-0.1))&(abs(yy)<(YYMax-0.1))
            l_x=int32(200*(xx+XXMax)/(2*XXMax));
            l_y=int32(200*(yy+YYMax)/(2*YYMax));
            Prob(l_x,l_y)=Prob(l_x,l_y)+1;
            countR=countR+1;
        end
        
        
    end
    %*****************************************************************************************************
    countR
    %?????? ???????????
    attrac = fopen('HistS20.dat','w');
    
    for ll_x=1:200
        for ll_y=1:200
            
            xx=(2*XXMax)*ll_x/200;
            xx=xx-XXMax;
            xx=xx/alfa;
            yy=(2*YYMax)*ll_y/200;
            yy=yy-YYMax;
            yy=yy/alfa;
            ppp=Prob(ll_x,ll_y)/countR;
            fprintf(attrac,'%8.4f %8.4f %8.4f\n',xx,yy,ppp);
        end
    end
    fclose(attrac);
    
    
end


%pcolor(l_x,l_y,Prob) % ?????? ????????????? ??????? P ?? ?????, ??????????? ????????? ??? ????????? l_x ? l_y.
%colorbar  % ??????? colorbar ????????? ? ???????? ??????? ???????????? ????? ???????.

% pr1=c_T(1,:).*(c_T(1,:)').';
% pr1_sr=pr1_sr+pr1;
