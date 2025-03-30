%% PSO-DBN
%% 清空环境变量
warning off             % 关闭警告信息
close all               % 关闭所有图窗
clear                   % 清空工作区变量
clc                     % 清空命令行

%% 添加路径
addpath('Toolbox\')      % 添加工具箱路径

%% 导入数据
P_train = xlsread('data','training set','B2:G181')';  % 导入训练集输入数据
T_train = xlsread('data','training set','H2:H181')';  % 导入训练集输出数据
P_test = xlsread('data','test set','B2:G45')';        % 导入测试集输入数据
T_test = xlsread('data','test set','H2:H45')';        % 导入测试集输出数据

%% 数据部分注意
M = size(P_train, 2);  % 训练样本数
N = size(P_test, 2);   % 测试样本数
outdim = 1;            % 输出维度

%% 数据归一化
[p_train, ps_input] = mapminmax(P_train, 0, 1);    % 输入数据归一化
p_test = mapminmax('apply', P_test, ps_input);     % 应用归一化参数到测试集输入数据
[t_train, ps_output] = mapminmax(T_train, 0, 1);   % 输出数据归一化
t_test = mapminmax('apply', T_test, ps_output);    % 应用归一化参数到测试集输出数据

%% 数据转置以适应模型
p_train = p_train'; 
p_test = p_test';
t_train = t_train'; 
t_test = t_test';

%% 初始化opts变量
opts.numepochs = 500;           % 训练次数
opts.batchsize = 12;            % 每次训练样本个数
opts.momentum  = 0;             % 动量参数
opts.alpha     = 0.05;          % 学习率

%% PSO优化的参数设置
dim = 4;                          % 优化参数的维数
lb  = [20, 20, 1000, 1.5];        % 优化参数下限
ub  = [80, 80, 3000, 3.0];        % 优化参数上限
numsum = length(ub) - 2;          % 隐藏层节点数
fun = @(x)fun(x, numsum, p_train, t_train, opts);  % 定义目标函数

pop = 6;                           % 粒子群规模
Max_iteration = 6;                 % 最大迭代次数
Vmax = 0.1;                        % 速度上限
Vmin = -0.1;                       % 速度下限

%% PSO算法优化DBN参数
tic;  % 添加计时开始
[Best_pos, Best_score, curve] = PSO(pop, Max_iteration, lb, ub, dim, fun, Vmax, Vmin); 
zbest = Best_pos;                 % 获取最佳位置

%% 打印最优参数
disp('------------------最优参数结果------------------')
disp(['隐藏层节点数: ', num2str(round(zbest(1:numsum)))])
disp(['学习率: ', num2str(zbest(end))])
disp(['迭代次数: ', num2str(round(zbest(numsum + 1)))])
disp(['最优适应度值: ', num2str(Best_score)])
disp('------------------------------------------------')

%% 提取最优参数
for i = 1:numsum
    dbn.sizes(i + 1) = round(zbest(i));  % 设置每层的节点数
end
dbn.sizes(1) = [];                       % 删除输入层节点数
lr = zbest(end);                         % 提取学习率
epochs = round(zbest(numsum + 1));       % 提取迭代次数

%% 训练DBN模型
dbn = dbnsetup(dbn, p_train, opts);      % 设置DBN模型
dbn = dbntrain(dbn, p_train, opts);      % 训练DBN模型

%% 添加输出层并反向微调
nn = dbnunfoldtonn(dbn, outdim);         % 展开DBN到神经网络
nn.activation_function = 'sigm';         % 设置激活函数
nn.learningRate = lr;                    % 设置学习率
nn.momentum = 0.5;                       % 设置动量参数

[nn, loss] = nntrain(nn, p_train, t_train, opts);  % 反向微调训练

%% 预测结果
t_sim1 = nnpredict(nn, p_train);         % 训练集预测
t_sim2 = nnpredict(nn, p_test);          % 测试集预测

%%  计算仿真时间并显示
simulationTime = toc;  % 计时结束
disp(['数值仿真总耗时: ', num2str(simulationTime), ' 秒']);

%% 数据反归一化
T_sim1 = mapminmax('reverse', t_sim1', ps_output);  % 训练集反归一化
T_sim2 = mapminmax('reverse', t_sim2', ps_output);  % 测试集反归一化

%% 绘制损失函数曲线
figure
plot(1:length(loss), loss, 'b-', 'LineWidth', 1)
xlim([1, length(loss)])
xlabel('迭代次数')
ylabel('误差损失')
legend('损失函数')
title('损失函数曲线')
grid on

%% 绘制收敛曲线
figure
plot(curve, 'linewidth', 1.5);
grid on
xlabel('迭代次数')
ylabel('适应度值')
title('PSO-DBN收敛曲线')

%% 测试集回归图与误差分析
figure;
plotregression(T_test, T_sim2, '回归图');  % 绘制回归图
figure;
ploterrhist(T_test - T_sim2, '误差直方图');  % 绘制误差直方图

%% 计算评估指标
error1 = sqrt(sum((T_sim1 - T_train).^2) ./ M);  % 训练集均方根误差RMSE
error2 = sqrt(sum((T_test - T_sim2).^2) ./ N);  % 测试集均方根误差RMSE
R1 = 1 - norm(T_train - T_sim1)^2 / norm(T_train - mean(T_train))^2;  % 训练集决定系数
R2 = 1 - norm(T_test - T_sim2)^2 / norm(T_test - mean(T_test))^2;     % 测试集决定系数
%均方误差 MSE
mse1 = sum((T_sim1 - T_train).^2)./M;
mse2 = sum((T_sim2 - T_test).^2)./N;

%RPD 剩余预测残差
SE1=std(T_sim1-T_train);
RPD1=std(T_train)/SE1;

SE=std(T_sim2-T_test);
RPD2=std(T_test)/SE;
%% 平均绝对误差MAE
MAE1 = mean(abs(T_train - T_sim1));
MAE2 = mean(abs(T_test - T_sim2));
%% 平均绝对百分比误差MAPE
MAPE1 = mean(abs((T_train - T_sim1)./T_train));
MAPE2 = mean(abs((T_test - T_sim2)./T_test));

%% 绘制预测结果对比图
figure
plot(1:M, T_train, 'r-*', 1:M, T_sim1, 'b-o', 'LineWidth', 1.5)
legend('真实值', 'PSO-DBN预测值')
xlabel('预测样本')
ylabel('预测结果')
title(['训练集预测结果对比 (R^2 =', num2str(R1), ', RMSE=', num2str(error1), ')'])

figure
plot(1:N, T_test, 'r-*', 1:N, T_sim2, 'b-o', 'LineWidth', 1.5)
legend('真实值', 'PSO-DBN预测值')
xlabel('预测样本')
ylabel('预测结果')
title(['测试集预测结果对比 (R^2 =', num2str(R2), ', RMSE=', num2str(error2), ')'])

%% 测试集误差图
figure  
ERROR3 = T_test - T_sim2;  % 计算误差
plot(ERROR3, 'b-*', 'LineWidth', 1.5)
xlabel('测试集样本编号')
ylabel('预测误差')
title('测试集预测误差')
grid on;
legend('PSO-DBN预测输出误差')

%% 打印评估指标
disp('-------------------测试集评估指标-------------------')
disp(['均方根误差RMSE: ', num2str(error2)])
disp(['决定系数R^2: ', num2str(R2)])


%% 绘制线性拟合图
%% 训练集拟合效果图
figure
plot(T_train,T_sim1,'*r');
xlabel('真实值')
ylabel('预测值')
string = {'训练集效果图';['R^2_c=' num2str(R1)  '  RMSEC=' num2str(error1) ]};
title(string)
hold on ;h=lsline;
set(h,'LineWidth',1,'LineStyle','-','Color',[1 0 1])
%% 预测集拟合效果图
figure
plot(T_test,T_sim2,'ob');
xlabel('真实值')
ylabel('预测值')
string1 = {'测试集效果图';['R^2_p=' num2str(R2)  '  RMSEP=' num2str(error2) ]};
title(string1)
hold on ;h=lsline();
set(h,'LineWidth',1,'LineStyle','-','Color',[1 0 1])
%% 求平均
R3=(R1+R2)./2;
error3=(error1+error2)./2;
%% 总数据线性预测拟合图
tsim=[T_sim1,T_sim2]';
S=[T_train,T_test]';
figure
plot(S,tsim,'ob');
xlabel('真实值')
ylabel('预测值')
string1 = {'所有样本拟合预测图';['R^2_p=' num2str(R3)  '  RMSEP=' num2str(error3) ]};
title(string1)
hold on ;h=lsline();
set(h,'LineWidth',1,'LineStyle','-','Color',[1 0 1])
%% 打印出评价指标
disp(['-----------------------误差计算--------------------------'])
disp(['训练集评价结果如下所示：'])
disp(['训练集平均绝对误差MAE为：',num2str(MAE1)])
disp(['训练集均方误差MSE为：       ',num2str(mse1)])
disp(['训练集均方根误差RMSE为：  ',num2str(error1)])
disp(['训练集决定系数R^2为：  ',num2str(R1)])
disp(['训练集剩余预测残差RPD为：  ',num2str(RPD1)])
disp(['训练集平均绝对百分比误差MAPE为：  ',num2str(MAPE1)])
disp(['测试集评价结果如下所示：'])
disp(['测试集平均绝对误差MAE为：',num2str(MAE2)])
disp(['测试集均方误差MSE为：       ',num2str(mse2)])
disp(['测试集均方根误差RMSE为：  ',num2str(error2)])
disp(['测试集决定系数R^2为：  ',num2str(R2)])
disp(['测试集剩余预测残差RPD为：  ',num2str(RPD2)])
disp(['测试集平均绝对百分比误差MAPE为：  ',num2str(MAPE2)])
grid

