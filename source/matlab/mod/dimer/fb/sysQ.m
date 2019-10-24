function dy = sysQ(t,x)
global N J E0 A0 w U g QEs QUs QJs Rs Ks phi0
T=2*pi/w; 

ft=sin(w*t+phi0);%t/T-floor(t/T);
ft=2*(ft>0)-1;

dy=(((E0+A0*sin(w*t+phi0))*QEs+U*QUs+J*QJs)+g/2*Rs)*x+g/2*Ks;
%dy=(((E0+A0*ft)*QEs+U*QUs+J*QJs)+g*Rs)*x+g*Ks;


end