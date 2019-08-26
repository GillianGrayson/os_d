function dy = sysQ(t,x)
global N J E0 A0 w U g QEs QUs QJs Rs Ks phi
   
dy=(((E0+1*A0*cos(w*t+phi))*QEs+(U+0*A0*cos(w*t+phi))*QUs+(J+0*A0*cos(w*t+phi))*QJs)+g*Rs)*x+g*Ks;

end

