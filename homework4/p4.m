clear; cvx_clear;

%% (a) & (b) & (c) 处理 v1
fprintf('--- 任务 (a, b, c): 处理 v1 ---\n');
v1 = [1.4; 0.5; -0.2; 2.0; 0.3];
n = length(v1);

% CVX 求解
cvx_begin quiet
    variable x_cvx(n)
    minimize( 0.5 * sum_square(x_cvx - v1) ) % 使用 sum_square 代替 norm(...)^2
    subject to
        x_cvx >= 0
        sum(x_cvx) == 1
cvx_end

% 算法求解
x_alg = project_simplex(v1);

% 报告结果
fprintf('CVX  解 (x_cvx): [%s]\n', num2str(x_cvx', '%.4f '));
fprintf('算法 解 (x_alg): [%s]\n', num2str(x_alg', '%.4f '));

% 计算指标
diff_norm = norm(x_cvx - x_alg, 2);
feas_res_sum = abs(sum(x_alg) - 1);
feas_res_min = min(x_alg);
% 驻留点残差：根据 KKT，x-v+mu=0，这里 mu 约等于 (v-x) 的均值
mu_est = mean(v1(x_alg > 0) - x_alg(x_alg > 0)); 
stat_res = norm(x_alg - v1 + mu_est, 2); 

fprintf('i.   二范数差异 ||x_cvx - x_alg||_2: %e\n', diff_norm);
fprintf('ii.  原可行性残差 |1^T x* - 1|: %e, min x_i: %f\n', feas_res_sum, feas_res_min);
fprintf('iii. 驻留点残差 (Stationarity): %e\n\n', stat_res);

%% (d) 处理 v2 (n=1000)
fprintf('--- 任务 (d): 大规模测试 (n=1000) ---\n');
rng(2026);
v2 = randn(1000, 1);

% 测试 CVX 时间
tic;
cvx_begin quiet
    variable x2_cvx(1000)
    minimize( 0.5 * sum_square(x2_cvx - v2) ) % 这里也改掉了！
    subject to
        sum(x2_cvx) == 1
        x2_cvx >= 0
cvx_end
t_cvx = toc;

% 测试排序算法时间
tic;
x2_alg = project_simplex(v2);
t_alg = toc;

fprintf('CVX 耗时: %.4f 秒\n', t_cvx);
fprintf('排序算法 耗时: %.4f 秒\n', t_alg);
fprintf('两者误差 ||x2_cvx - x2_alg||: %e\n', norm(x2_cvx - x2_alg));

%% --- 函数定义部分必须放在文件末尾 ---
function x = project_simplex(v)
    n = length(v);
    u = sort(v, 'descend');
    sv = cumsum(u);
    % 寻找满足 u_k - (1/k)*(sum(u_1...u_k) - 1) > 0 的最大 k
    rho = find(u - (sv - 1) ./ (1:n)' > 0, 1, 'last');
    theta = (sv(rho) - 1) / rho;
    x = max(v - theta, 0);
end