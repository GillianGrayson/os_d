clear all;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/purity_final.txt', data_path);
pur_fin_data = importdata(fn);
pur_fin = pur_fin_data(1)

fn = sprintf('%s/purity_avg.txt', data_path);
pur_avg_data = importdata(fn);
pur_avg = pur_avg_data(1)
