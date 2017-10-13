clear all;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/eigs.txt', data_path);
data = importdata(fn);

fig = figure;
propertyeditor(fig);

scatter(data(:, 1), data(:, 2), 'LineWidth', 2, 'MarkerSize', 8);

hold all;

