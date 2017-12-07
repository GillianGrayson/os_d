clear all;

N = 10;
J = -1;

for U = 0.01 : 0.01 : 1.00
    
    U = U
    
    g = 0.1;
    
    to_data_path = '../data/qj_input';
    
    %????????? ??? E(t)
    e0 = 1.0;
    e1 = -1.0;
    tau1 = pi;
    T = 2*pi;
    
    E = 0.0;
    amplitude = 0.0;
    period = T;
    
    out_suffix = sprintf('E%0.4f_T%0.4f_A%0.4f_N%d_U%0.4f_J%0.4f_g%0.4f.bin', ...
        E, ...
        period, ...
        amplitude, ...
        N, ...
        U, ...
        J, ...
        g);
    
    %?????? ?????? a ? a_
    U = U / (N-1);
    g = g / (N-1);
    
    %????????????? ?????????
    k = 1;
    A1 = zeros(N);
    A2 = zeros(N);
    for i=1:N
        A1(i,i) = i - 1;
        A2(i,i) = N - i;
    end
    
    %??????????? ???????????? H(t)=H1+F(t)*H2
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
    
    hamiltonian = H1 + E * H2;
    
    H1 =  H1 - (i/2) * (g * V' * V);
    
    E0 = E + amplitude * e0;
    E1 = E + amplitude * e1;
    t1 = tau1; %????? ??????????, ????? F ?? ????? 0
    t2 = T - t1; %????? ??????????, ????? F ????? 0
    
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
    split = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
       
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
            ch1;
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
    split = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
    
    
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
            ch2;
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
    
    file_name_data = sprintf('%s/aux_data_%s', ...
        to_data_path , ...
        out_suffix);
    fid = fopen(file_name_data,'wb');
    
    res_H = zeros(2*N*N,1);
    cur_id = 1;
    for state_id_1 = 1:N
        for state_id_2 = 1:N
            res_H(cur_id) = real(hamiltonian(state_id_1, state_id_2));
            res_H(cur_id+1) = imag(hamiltonian(state_id_1, state_id_2));
            cur_id = cur_id + 2;
        end
    end
    fwrite(fid, res_H, 'double');
    res_H = 0;
    fclose(fid);
    
end