function dy = sysLind(t,x)
global N A0 w A phi0 Ps Pd
T=2*pi/w; 


ft=sin(w*t+phi0);%t/T-floor(t/T);
%ft=2*(ft>0)-1;

dy=(Ps+Pd*A0*ft)*x;

%for i=1:N+1
%    Rho(:,i)=x(1+(i-1)*(N+1):i*(N+1),1);
%end


%H=((E0+1*A0*ft)*HE+U/N*HU+J*HJ);
%dy2=-sqrt(-1)*(H*Rho-Rho*H)+g/2*(2*A*Rho*A'-Rho*A'*A-A'*A*Rho);


%   for i=1:N+1
%    dy(1+(i-1)*(N+1):i*(N+1),1)=dy2(:,i);
%   end

   
end

