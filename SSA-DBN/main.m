%% SSA-DBN
%%  清空环境变量
warning off             % 关闭报警信息
close all               % 关闭开启的图窗
clear                   % 清空变量
clc                     % 清空命令行

%%  添加路径
addpath('Toolbox\')

%%  导入数据
P_train = xlsread('data','training set','B2：G181')';
T_train= xlsread('data','training set','H2：H181')';
% 测试集——44个样本
P_test=xlsread('data','test set','B2:G45')';
T_test=xlsread('data','test set','H2:H45')';

%% 数据部分注意
% （M / batchsize = 整数）  opts.batchsize= 12;   
% 输入特征维度
%% 数据部分注意

M = size(P_train, 2);
N = size(P_test, 2);
outdim = 1;                                  % 最后一列为输出
%%  数据归一化
[p_train, ps_input] = mapminmax(P_train, 0, 1);
p_test = mapminmax('apply', P_test, ps_input);

[t_train, ps_output] = mapminmax(T_train, 0, 1);
t_test = mapminmax('apply', T_test, ps_output);

%%  转置以适应模型
p_train = p_train'; p_test = p_test';
t_train = t_train'; t_test = t_test';

%%  预训练的参数设置
opts.numepochs = 500;                   % 训练次数
opts.batchsize = 12;                    % 每次训练样本个数 需满足：（M / batchsize = 整数）
opts.momentum  = 0;                     % 学习率的动量
opts.alpha     = 0.05;                  % 学习率

%% 参数设置
dim = 4;                              % 优化参数个数
lb  = [20,  20, 1000,  1.5];       % 优化参数目标下限
ub  = [80,  80, 3000,  3.0];       % 优化参数目标上限
% ------ 前 numsum 个数为隐藏层节点边界和速度，最后一个为反向学习率，倒数第二个为迭代次数 -------
numsum = length(ub) - 2;              % 隐藏层层数
fun = @(x)fun(x,numsum, p_train, t_train, opts);                   % 目标函数
pop = 6;                              % 麻雀数量
Max_iteration = 6;                   % 最大迭代次数   

%%  优化算法
tic;  % 添加计时开始
[Best_score,Best_pos, curve] = SSA(pop, Max_iteration, lb, ub, dim, fun); 
zbest = Best_pos;
%%  提取最优参数
for i = 1 : numsum
    dbn.sizes(i + 1) = round(zbest(i));
end
dbn.sizes(1) = [];                      %
lr = zbest(end);                        % 学习率
epochs = round(zbest(numsum + 1));      % 迭代次数

%%  训练模型
dbn = dbnsetup(dbn, p_train, opts);     % 建立模型
dbn = dbntrain(dbn, p_train, opts);     % 训练模型

%%  训练权重移植，添加输出层
nn = dbnunfoldtonn(dbn, outdim);

%%  反向调整网络
% ----- 参数修改时，请查看fun函数中的参数设置，并作出对应修改 ---------
opts.numepochs          = epochs;       % 反向微调次数
opts.batchsize          = 12;           % 每次反向微调样本数 需满足：（M / batchsize = 整数）

nn.activation_function  = 'sigm';       % 激活函数
nn.learningRate         = lr;           % 学习率
nn.momentum             = 0.5;          % 动量参数
nn.scaling_learningRate = 1;            % 学习率的比例因子

[nn, loss] = nntrain(nn, p_train, t_train, opts);  % 训练

%%  预测 
t_sim1 = nnpredict(nn, p_train); 
t_sim2 = nnpredict(nn, p_test );

%%  计算仿真时间并显示
simulationTime = toc;  % 计时结束
disp(['数值仿真总耗时: ', num2str(simulationTime), ' 秒']);
%%  数据反归一化
T_sim1 = mapminmax('reverse', t_sim1', ps_output);
T_sim2 = mapminmax('reverse', t_sim2', ps_output);

%%  绘制损失函数曲线
figure
plot(1: length(loss), loss, 'b-', 'LineWidth', 1)
xlim([1, length(loss)])
xlabel('迭代次数')
ylabel('误差损失')
legend('损失函数')
title('损失函数')
grid

%% 绘图并保存最佳变量
figure
plot( curve,'linewidth',1.5);
grid on
xlabel('迭代次数')
ylabel('适应度函数')
title('SSA-DBN收敛曲线')
%% 测试集结果
figure;
plotregression(T_test,T_sim2,['回归图']);
figure;
ploterrhist(T_test-T_sim2,['误差直方图']);
%%  均方根误差 RMSE
error1 = sqrt(sum((T_sim1 - T_train).^2)./M);
error2 = sqrt(sum((T_test - T_sim2).^2)./N);

%%
%决定系数
R1 = 1 - norm(T_train - T_sim1)^2 / norm(T_train - mean(T_train))^2;
R2 = 1 - norm(T_test -  T_sim2)^2 / norm(T_test -  mean(T_test ))^2;

%%
%均方误差 MSE
mse1 = sum((T_sim1 - T_train).^2)./M;
mse2 = sum((T_sim2 - T_test).^2)./N;
%%
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
%%  训练集绘图
figure
%plot(1:M,T_train,'r-*',1:M,T_sim1,'b-o','LineWidth',1)
plot(1:M,T_train,'r-*',1:M,T_sim1,'b-o','LineWidth',1.5)
legend('真实值','SSA-DBN预测值')
xlabel('预测样本')
ylabel('预测结果')
string={'训练集预测结果对比';['(R^2 =' num2str(R1) ' RMSE= ' num2str(error1) ' MSE= ' num2str(mse1) ' RPD= ' num2str(RPD1) ')' ]};
title(string)
%% 预测集绘图
figure
plot(1:N,T_test,'r-*',1:N,T_sim2,'b-o','LineWidth',1.5)
legend('真实值','SSA-DBN预测值')
xlabel('预测样本')
ylabel('预测结果')
string={'测试集预测结果对比';['(R^2 =' num2str(R2) ' RMSE= ' num2str(error2)  ' MSE= ' num2str(mse2) ' RPD= ' num2str(RPD2) ')']};
title(string)

%% 测试集误差图
figure  
ERROR3=T_test-T_sim2;
plot(T_test-T_sim2,'b-*','LineWidth',1.5)
xlabel('测试集样本编号')
ylabel('预测误差')
title('测试集预测误差')
grid on;
legend('SSA-DBN预测输出误差')
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