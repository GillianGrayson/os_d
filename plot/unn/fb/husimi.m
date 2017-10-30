function y = husimi(theta, phi, rho)

sizes = size(rho);
N = sizes(1)-1;
Ngrid = length(theta);
Mgrid = length(phi);
cos_theta = cos(theta/2);
sin_theta = sin(theta/2);
exp_iphi = exp(sqrt(-1)*phi);

y = zeros(Ngrid,Mgrid);

for i=1:Ngrid
    for j=1:Mgrid
        i = i
        j = j
        
        for n=0:N
            coherent(n+1,1) = sqrt(nchoosek(N,n))*cos_theta(i)^n*(sin_theta(i)*exp_iphi(j))^(N-n);
        end
        y(i,j) = (coherent')*rho*coherent;
    end
end

y=y/sum(sum(y));

end