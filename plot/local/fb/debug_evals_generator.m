clear all;

N = 20;
NS = (N+1) * (N+1) - 1;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/eigs.txt', data_path);
curr_data = importdata(fn);

scatter(curr_data(:, 1), curr_data(:, 2))
    
max(curr_data(:, 1))
