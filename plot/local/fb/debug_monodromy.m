clear all;

N = 5;
NS = (N+1) * (N+1) - 1;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';


mndrm_init = zeros(NS, NS);
mndrm_fin = zeros(NS, NS);

for mult_id = 1 : NS
    
    fn = sprintf('%s/vec_%d_init.txt', data_path, mult_id - 1);
    curr_data = importdata(fn);
    
    for st_id = 1 : NS
        mndrm_init(st_id, mult_id) = curr_data(st_id, 1);
    end
    
    fn = sprintf('%s/vec_%d_fin.txt', data_path, mult_id - 1);
    curr_data = importdata(fn);
    
    for st_id = 1 : NS
        mndrm_fin(st_id, mult_id) = curr_data(st_id, 1);
    end

end

max(max(mndrm_init - eye(NS)))

evals = eig(mndrm_fin);