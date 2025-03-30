%% DBN
%% 清空环境变量
warning off             % 关闭报警信息
close all               % 关闭开启的图窗
clear                   % 清空变量
clc                     % 清空命令行

%% 导入数据
P_train = xlsread('data', 'training set', 'B2:G181')';  % 导入训练集输入数据
T_train = xlsread('data', 'training set', 'H2:H181')';  % 导入训练集输出数据
P_test  = xlsread('data', 'test set', 'B2:G45')';       % 导入测试集输入数据
T_test  = xlsread('data', 'test set', 'H2:H45')';       % 导入测试集输出数据

M = size(P_train, 2);   % 训练集样本个数
N = size(P_test, 2);    % 测试集样本个数

%% 数据归一化
[p_train, ps_input] = mapminmax(P_train, 0, 1);
p_test = mapminmax('apply', P_test, ps_input);

[t_train, ps_output] = mapminmax(T_train, 0, 1);
t_test = mapminmax('apply', T_test, ps_output);

%% 转置以适应模型
p_train = p_train'; p_test = p_test';
t_train = t_train'; t_test = t_test';

%% 模型预训练
dbn.sizes      = [50, 50];           % 隐藏层节点数
opts.numepochs = 500;                % 训练次数
opts.batchsize = 12;                 % 每次训练样本个数，需满足：（M / batchsize = 整数）
opts.momentum  = 0;                  % 动量参数
opts.alpha     = 0.05;                % 学习率

dbn = dbnsetup(dbn, p_train, opts);  % 建立模型
dbn = dbntrain(dbn, p_train, opts);  % 训练模型

%% 训练权重移植，添加输出层
tic;  % 添加计时开始
nn = dbnunfoldtonn(dbn, size(T_train, 1));

%% 反向调整网络
opts.numepochs = 2000;               % 反向微调次数
opts.batchsize = 12;                 % 每次反向微调样本数，需满足：（M / batchsize = 整数）

nn.activation_function = 'sigm';     % 激活函数
nn.learningRate = 0.4;                 % 学习率
nn.momentum = 0.5;                   % 动量参数
nn.scaling_learningRate = 1;         % 学习率的比例因子

[nn, loss] = nntrain(nn, p_train, t_train, opts); % 反向微调训练

%% 模型预测
t_sim1 = nnpredict(nn, p_train);
t_sim2 = nnpredict(nn, p_test);

%%  计算仿真时间并显示
simulationTime = toc;  % 计时结束
disp(['数值仿真总耗时: ', num2str(simulationTime), ' 秒']);

%% 数据反归一化
T_sim1 = mapminmax('reverse', t_sim1, ps_output);
T_sim2 = mapminmax('reverse', t_sim2, ps_output);

%% 测试集结果回归图与误差分析
figure;
plotregression(T_test, T_sim2, '回归图');
figure;
ploterrhist(T_test - T_sim2, '误差直方图');

%% 均方根误差
error1 = sqrt(sum((T_sim1' - T_train).^2) ./ M);
error2 = sqrt(sum((T_sim2' - T_test ).^2) ./ N);

%% 绘制损失函数曲线
figure
plot(1:length(loss), loss, 'b-', 'LineWidth', 1)
xlim([1, length(loss)])
xlabel('迭代次数')
ylabel('误差损失')
legend('损失函数')
title('损失函数')
grid

%% 绘图
figure
plot(1:M, T_train, '-s', 1:M, T_sim1, '-o', 'LineWidth', 1)
legend('真实值', '预测值')
xlabel('预测样本')
ylabel('预测结果')
string = {'训练集预测结果对比'; ['RMSE=' num2str(error1)]};
title(string)
xlim([1, M])
grid

figure
plot(1:N, T_test, '-s', 1:N, T_sim2, '-o', 'LineWidth', 1)
legend('真实值', '预测值')
xlabel('预测样本')
ylabel('预测结果')
string = {'测试集预测结果对比'; ['RMSE=' num2str(error2)]};
title(string)
xlim([1, N])
grid

%% 计算评估指标
% R2
R1 = 1 - norm(T_train - T_sim1')^2 / norm(T_train - mean(T_train))^2;
R2 = 1 - norm(T_test - T_sim2')^2 / norm(T_test - mean(T_test))^2;
% MAE
mae1 = sum(abs(T_sim1' - T_train)) ./ M;
mae2 = sum(abs(T_sim2' - T_test)) ./ N;
% MSE
mse1 = sum((T_sim1' - T_train).^2) ./ M;
mse2 = sum((T_sim2' - T_test).^2) ./ N;
% MAPE
mape1 = sum(abs((T_sim1' - T_train) ./ T_train)) ./ M * 100;
mape2 = sum(abs((T_sim2' - T_test) ./ T_test)) ./ N * 100;
% RMSE
rmse1 = sqrt(mse1);
rmse2 = sqrt(mse2);
% RPD（剩余预测残差）
rpd1 = std(T_train) / std(T_train - T_sim1');
rpd2 = std(T_test) / std(T_test - T_sim2');

%% 绘制线性拟合图
% 训练集拟合效果图
figure
plot(T_train, T_sim1, '*r');
xlabel('真实值');
ylabel('预测值');
string = {'训练集效果图'; ['R^2_c=' num2str(R1) '  RMSEC=' num2str(error1)]};
title(string);
hold on;
h = lsline; % 添加线性拟合线
set(h, 'LineWidth', 1, 'LineStyle', '-', 'Color', [1 0 1]); % 设置拟合线样式

% 测试集拟合效果图
figure
plot(T_test, T_sim2, 'ob');
xlabel('真实值');
ylabel('预测值');
string1 = {'测试集效果图'; ['R^2_p=' num2str(R2) '  RMSEP=' num2str(error2)]};
title(string1);
hold on;
h = lsline(); % 添加线性拟合线
set(h, 'LineWidth', 1, 'LineStyle', '-', 'Color', [1 0 1]); % 设置拟合线样式


%% 打印评估指标
disp('-------------------评估指标-------------------')
disp(['训练集 R^2: ', num2str(R1)])
disp(['测试集 R^2: ', num2str(R2)])
disp(['训练集 MSE: ', num2str(mse1)])
disp(['测试集 MSE: ', num2str(mse2)])
disp(['训练集 MAE: ', num2str(mae1)])
disp(['测试集 MAE: ', num2str(mae2)])
disp(['训练集 MAPE: ', num2str(mape1), '%'])
disp(['测试集 MAPE: ', num2str(mape2), '%'])
disp(['训练集 RMSE: ', num2str(rmse1)])
disp(['测试集 RMSE: ', num2str(rmse2)])
disp(['训练集 RPD: ', num2str(rpd1)])
disp(['测试集 RPD: ', num2str(rpd2)])
disp('---------------------------------------------')


