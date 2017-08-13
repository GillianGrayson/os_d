N = 100;
J = 1;
for U = 0.5:0.1:0.5
    U = U
    
    g = 0.1;
    
    to_data_path = '../data/qj_input';
    
    out_suffix = sprintf('N%d_U%0.4f.bin', ...
        N, ...
        U);
    
    %параметры для E(t)
    e0 = 1.0;
    e1 = -1.0;
    tau1 = pi;
    T = 2*pi;
    
    %размер матриц a и a_
    U=U/(N-1);
    g=g/(N-1);
    
    %диссипативные операторы
    k = 1;
    A1 = zeros(N);
    A2 = zeros(N);
    for i=1:N
        A1(i,i) = i - 1;
        A2(i,i) = N - i;
    end
    
    %эффективный гамильтониан H(t)=H1+F(t)*H2
    H1 = zeros(N);
    H2 = A2 - A1;
    
    %H=-J(b1'*b2+b2'*b1)+U/2*(n1(n1-1)+n2*(n2-1))(+-)e0*(n2-n1)
    
    H1 = U*2*(A1*(A1-eye(N)) + A2*(A2-eye(N)));
    
    for i=1:N-1
        %     if (i+1<=Nn)
        %         H1(i, i+1) = H1(i, i+1)-J*sqrt((Nn-i)*(i+1)); %b1'*b2
        %     end
        %     if (i-1>0)
        %         H1(i, i-1) = H1(i, i-1)-J*sqrt((i)*(Nn-i+1)); %b2'*b1
        %     end
        
        H1(i+1,i)=H1(i+1,i)-J*sqrt((N-i)*(i));
        H1(i,i+1)=H1(i,i+1)-J*sqrt((i)*(N-i));
    end
    
    %V=(b1'+b2')(b1-b2)=b1'*b1-b2'*b2-b1'b2+b2'*b1=n1-n2+b2'*b1-b1'*b2
    V = A1 - A2;
    for i=1:N-1
        %     if (i+1<=N)
        %         V(i, i+1) = V(i, i+1)-sqrt((N-i+1)*(i)); %b1'*b2
        %     end
        %     if (i-1>0)
        %         V(i, i-1) = V(i, i-1)+sqrt((i)*(N-i+1)); %b2'*b1
        %     end
        %
        V(i+1,i)=V(i+1,i)-sqrt((N-i)*(i));
        V(i,i+1)=V(i,i+1)+sqrt((i)*(N-i));
    end
    
    
    i = sqrt(-1);
    
    H1 =  H1 - (i/2) * (g * V' * V);
    
    E0 = 1 + 1.5 * e0;
    E1 = 1 + 1.5 * e1;
    t1 = tau1; %длина промежутка, когда F не равно 0
    t2 = T- t1; %длина промежутка, когда F равно 0
    
    res = zeros(2*N*N,1);
    
    file_name_data = sprintf('%s/main_data_%s', ...
        to_data_path , ...
        out_suffix);
    fid = fopen(file_name_data,'wb');
    
    k=2;
    fwrite(fid, T, 'double');
    fwrite(fid, N, 'int');
    fwrite(fid, k, 'int');
    
    dt = t1;
    deep = 16;
    split = [64, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
    
    
    
    fwrite(fid, deep, 'int');
    for j = 1:deep
        dt = dt/split(j);
        
        fwrite(fid, dt, 'double');
        fwrite(fid, split(j), 'int');
        G1 = expm(-i*(H1+E0*H2)*dt);
        
        H3 = (-i*(H1+E0*H2)*dt);
        [R, D] = eig(G1);
        ph = ones (N, 1);
        nr = norm(ph);
        ph1 = G1*ph;
        ch1 = norm(ph1)/norm(ph);
        if (j==1)
            ch1
        end
        
        m=1;
        for l = 1:N
            for k=1:N
                res(m)= real(G1(l,k));
                res(m+1)=imag(G1(l,k));
                m=m+2;
            end
        end
        fwrite(fid, res, 'double');
    end
    %
    dt = t2;
    deep = 16;
    split = [64, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
    
    
    fwrite(fid, deep, 'int');
    for j = 1:deep
        dt = dt/split(j);
        
        fwrite(fid, dt, 'double');
        fwrite(fid, split(j), 'int');
        
        G2=expm(-i*(H1+E1*H2)*dt);
        
        H3 = (-i*(H1+E1*H2)*dt);
        [R, D] = eig(G2);
        ph = ones (N, 1);
        nr = norm(ph);
        ph1 = G2*ph;
        ch2= norm(ph1)/norm(ph);
        if (j==1)
            ch2
        end
        
        m=1;
        for l = 1:N
            for k=1:N
                res(m)= real(G2(l,k));
                res(m+1)=imag(G2(l,k));
                m=m+2;
            end
        end
        fwrite(fid, res, 'double');
    end
    %
    
    k = 1;
    fwrite(fid, k, 'int');
    
    One=1.0;
    
    fwrite(fid, One, 'double');
    m=1;
    for l = 1:N
        for k=1:N
            res(m)= real(V(l,k));
            res(m+1)=imag(V(l,k));
            m=m+2;
        end
    end
    fwrite(fid, res, 'double');
    
    fclose(fid);
    
end