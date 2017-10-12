clear all;

data_path = '../../../source/cpp/CQdiss/CQdiss_fbasis';

fn = sprintf('%s/multiplicators.txt', data_path);
data = importdata(fn);

fig = figure;
propertyeditor(fig);

circle(0, 0, 1)
hold all;
xlim([-2 2])
ylim([-2 2])

scatter(data(:, 1), data(:, 2), 'LineWidth', 2, 'MarkerSize', 8);

hold all;




function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end