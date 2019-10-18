clearvars h

num_points = 500;

h = animatedline('Color', 'm', 'MaximumNumPoints', num_points, 'LineWidth', 1);


for d_id = 1:global_size
    curr_re = all_re(d_id, 1);
    curr_im = all_im(d_id, 1);

    addpoints(h, curr_re, curr_im);
    drawnow
end




