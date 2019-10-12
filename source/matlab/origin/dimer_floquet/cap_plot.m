detltas = zeros(size(Us, 2), 1);

for U_id = 1:size(Us, 2)
    eigv = evals(:, U_id);
    abs_evals = sort(abs(eigv));
    
    detltas(U_id) = 1 - abs_evals(end - 1);
end

a = 0
hLine = plot(Us, detltas, 'LineWidth', 2);
set(gca, 'FontSize', 30);
xlabel('$U$', 'Interpreter', 'latex');
set(gca, 'FontSize', 30);
ylabel('$\Delta$', 'Interpreter', 'latex');
hold all;


% figure;
% hLine = plot(real(eigvals(:, U_id)), imag(eigvals(:, U_id)), '.', 'LineWidth', 2, 'MarkerSize', 10);
% set(gca, 'FontSize', 30);
% xlabel('$re$', 'Interpreter', 'latex');
% set(gca, 'FontSize', 30);
% ylabel('$im$', 'Interpreter', 'latex');
% xlim([-1.1 1.1])
% ylim([-1.1 1.1])
% hold all;
% ang=0:0.01:2*pi; 
% xp=cos(ang);
% yp=sin(ang);
% plot(xp,yp);
