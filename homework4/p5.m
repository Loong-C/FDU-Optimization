clear; clc;

% 数据定义
t = [0; 1; 2; 3; 4; 5];
c = [1; 2; 3; 4; 5; 9];
A = [ones(6,1), t];

%% (3) 求解 L1 范数拟合
cvx_begin quiet
    variables x1(2)
    minimize( norm(A*x1 - c, 1) )
cvx_end

fprintf('L1 拟合结果 (a*, b*): [%.4f, %.4f]\n', x1(1), x1(2));
fprintf('L1 最优值: %.4f\n', cvx_optval);

%% (4) 求解 L2 最小二乘拟合
cvx_begin quiet
    variables x2(2)
    minimize( norm(A*x2 - c, 2) ) % 或者用 sum_square(A*x2 - c)
cvx_end

fprintf('L2 拟合结果 (a*, b*): [%.4f, %.4f]\n', x2(1), x2(2));
fprintf('L2 最优值 (范数): %.4f\n', cvx_optval);

% 绘图
figure;
plot(t, c, 'ko', 'MarkerFaceColor', 'k'); hold on;
t_plot = linspace(0, 5, 100);
plot(t_plot, x1(1) + x1(2)*t_plot, 'r-', 'LineWidth', 1.5); % L1 线
plot(t_plot, x2(1) + x2(2)*t_plot, 'b--', 'LineWidth', 1.5); % L2 线
legend('原始数据', 'L1 拟合 (Robust)', 'L2 拟合 (Least-squares)');
grid on;
title('L1 vs L2 拟合对比');