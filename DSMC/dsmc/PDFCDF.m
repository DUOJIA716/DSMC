% MATLAB代码：线性与指数分布的PDF和CDF对比

% 参数设置
alpha_max = 0.63; % 最大填充率 alpha_max
lambda = 1;       % 指数分布的衰减参数 lambda

% 定义 alpha 的范围
alpha = linspace(0, alpha_max, 1000); % 从0到alpha_max等间隔取1000个点

% 线性分布的PDF和CDF
% 线性PDF
p_alpha_linear = (2 / alpha_max) * (1 - alpha / alpha_max);
% 线性CDF
F_alpha_linear = (2 / alpha_max) * (alpha - (alpha.^2) / (2 * alpha_max));

% 指数分布的PDF和CDF
% 计算归一化常数A
numerator_A = 1; % A的分子
denominator_A = (1 - exp(-lambda * alpha_max)) * (1 - 1 / alpha_max) + exp(-lambda * alpha_max); % A的分母
A = numerator_A / denominator_A; % 归一化常数A
% 指数PDF
p_alpha_exp = A .* (1 - alpha / alpha_max) .* exp(-lambda * alpha);
% 指数CDF
F_alpha_exp = A .* ((1 - exp(-lambda * alpha)) - (1 / alpha_max) * (-alpha .* exp(-lambda * alpha) + (1 - exp(-lambda * alpha))));

% 绘制PDF对比图
figure;
plot(alpha, p_alpha_linear, 'b--', 'LineWidth', 2); hold on; % 绘制线性PDF
plot(alpha, p_alpha_exp, 'r', 'LineWidth', 2);              % 绘制指数PDF
xlabel('\alpha (填充率)'); % 横轴标签
ylabel('PDF p(\alpha)');  % 纵轴标签
title('线性与指数分布的PDF对比'); % 标题
legend('线性PDF', '指数PDF'); % 图例
grid on; % 添加网格

% 绘制CDF对比图
figure;
plot(alpha, F_alpha_linear, 'b--', 'LineWidth', 2); hold on; % 绘制线性CDF
plot(alpha, F_alpha_exp, 'r', 'LineWidth', 2);              % 绘制指数CDF
xlabel('\alpha (填充率)'); % 横轴标签
ylabel('CDF F(\alpha)');  % 纵轴标签
title('线性与指数分布的CDF对比'); % 标题
legend('线性CDF', '指数CDF'); % 图例
grid on; % 添加网格

% 显示指数分布的归一化常数A
disp('指数分布的归一化常数 A:');
disp(A);
