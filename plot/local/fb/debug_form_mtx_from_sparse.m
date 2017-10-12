N=21;

curr_rho = zeros(N, N);
for d_id = 1:size(rho_data, 1)
    curr_row = rho_data(d_id, 1);
    curr_col = rho_data(d_id, 2);
    curr_rho(curr_row, curr_col) = rho_data(d_id, 3) + sqrt(-1) * rho_data(d_id, 4);
end
