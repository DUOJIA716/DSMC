%% SSA-DBN
%%  ��ջ�������
warning off             % �رձ�����Ϣ
close all               % �رտ�����ͼ��
clear                   % ��ձ���
clc                     % ���������

%%  ���·��
addpath('Toolbox\')

%%  ��������
P_train = xlsread('data','training set','B2��G181')';
T_train= xlsread('data','training set','H2��H181')';
% ���Լ�����44������
P_test=xlsread('data','test set','B2:G45')';
T_test=xlsread('data','test set','H2:H45')';

%% ���ݲ���ע��
% ��M / batchsize = ������  opts.batchsize= 12;   
% ��������ά��
%% ���ݲ���ע��

M = size(P_train, 2);
N = size(P_test, 2);
outdim = 1;                                  % ���һ��Ϊ���
%%  ���ݹ�һ��
[p_train, ps_input] = mapminmax(P_train, 0, 1);
p_test = mapminmax('apply', P_test, ps_input);

[t_train, ps_output] = mapminmax(T_train, 0, 1);
t_test = mapminmax('apply', T_test, ps_output);

%%  ת������Ӧģ��
p_train = p_train'; p_test = p_test';
t_train = t_train'; t_test = t_test';

%%  Ԥѵ���Ĳ�������
opts.numepochs = 500;                   % ѵ������
opts.batchsize = 12;                    % ÿ��ѵ���������� �����㣺��M / batchsize = ������
opts.momentum  = 0;                     % ѧϰ�ʵĶ���
opts.alpha     = 0.05;                  % ѧϰ��

%% ��������
dim = 4;                              % �Ż���������
lb  = [20,  20, 1000,  1.5];       % �Ż�����Ŀ������
ub  = [80,  80, 3000,  3.0];       % �Ż�����Ŀ������
% ------ ǰ numsum ����Ϊ���ز�ڵ�߽���ٶȣ����һ��Ϊ����ѧϰ�ʣ������ڶ���Ϊ�������� -------
numsum = length(ub) - 2;              % ���ز����
fun = @(x)fun(x,numsum, p_train, t_train, opts);                   % Ŀ�꺯��
pop = 6;                              % ��ȸ����
Max_iteration = 6;                   % ����������   

%%  �Ż��㷨
tic;  % ��Ӽ�ʱ��ʼ
[Best_score,Best_pos, curve] = SSA(pop, Max_iteration, lb, ub, dim, fun); 
zbest = Best_pos;
%%  ��ȡ���Ų���
for i = 1 : numsum
    dbn.sizes(i + 1) = round(zbest(i));
end
dbn.sizes(1) = [];                      %
lr = zbest(end);                        % ѧϰ��
epochs = round(zbest(numsum + 1));      % ��������

%%  ѵ��ģ��
dbn = dbnsetup(dbn, p_train, opts);     % ����ģ��
dbn = dbntrain(dbn, p_train, opts);     % ѵ��ģ��

%%  ѵ��Ȩ����ֲ����������
nn = dbnunfoldtonn(dbn, outdim);

%%  �����������
% ----- �����޸�ʱ����鿴fun�����еĲ������ã���������Ӧ�޸� ---------
opts.numepochs          = epochs;       % ����΢������
opts.batchsize          = 12;           % ÿ�η���΢�������� �����㣺��M / batchsize = ������

nn.activation_function  = 'sigm';       % �����
nn.learningRate         = lr;           % ѧϰ��
nn.momentum             = 0.5;          % ��������
nn.scaling_learningRate = 1;            % ѧϰ�ʵı�������

[nn, loss] = nntrain(nn, p_train, t_train, opts);  % ѵ��

%%  Ԥ�� 
t_sim1 = nnpredict(nn, p_train); 
t_sim2 = nnpredict(nn, p_test );

%%  �������ʱ�䲢��ʾ
simulationTime = toc;  % ��ʱ����
disp(['��ֵ�����ܺ�ʱ: ', num2str(simulationTime), ' ��']);
%%  ���ݷ���һ��
T_sim1 = mapminmax('reverse', t_sim1', ps_output);
T_sim2 = mapminmax('reverse', t_sim2', ps_output);

%%  ������ʧ��������
figure
plot(1: length(loss), loss, 'b-', 'LineWidth', 1)
xlim([1, length(loss)])
xlabel('��������')
ylabel('�����ʧ')
legend('��ʧ����')
title('��ʧ����')
grid

%% ��ͼ��������ѱ���
figure
plot( curve,'linewidth',1.5);
grid on
xlabel('��������')
ylabel('��Ӧ�Ⱥ���')
title('SSA-DBN��������')
%% ���Լ����
figure;
plotregression(T_test,T_sim2,['�ع�ͼ']);
figure;
ploterrhist(T_test-T_sim2,['���ֱ��ͼ']);
%%  ��������� RMSE
error1 = sqrt(sum((T_sim1 - T_train).^2)./M);
error2 = sqrt(sum((T_test - T_sim2).^2)./N);

%%
%����ϵ��
R1 = 1 - norm(T_train - T_sim1)^2 / norm(T_train - mean(T_train))^2;
R2 = 1 - norm(T_test -  T_sim2)^2 / norm(T_test -  mean(T_test ))^2;

%%
%������� MSE
mse1 = sum((T_sim1 - T_train).^2)./M;
mse2 = sum((T_sim2 - T_test).^2)./N;
%%
%RPD ʣ��Ԥ��в�
SE1=std(T_sim1-T_train);
RPD1=std(T_train)/SE1;

SE=std(T_sim2-T_test);
RPD2=std(T_test)/SE;
%% ƽ���������MAE
MAE1 = mean(abs(T_train - T_sim1));
MAE2 = mean(abs(T_test - T_sim2));
%% ƽ�����԰ٷֱ����MAPE
MAPE1 = mean(abs((T_train - T_sim1)./T_train));
MAPE2 = mean(abs((T_test - T_sim2)./T_test));
%%  ѵ������ͼ
figure
%plot(1:M,T_train,'r-*',1:M,T_sim1,'b-o','LineWidth',1)
plot(1:M,T_train,'r-*',1:M,T_sim1,'b-o','LineWidth',1.5)
legend('��ʵֵ','SSA-DBNԤ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ����')
string={'ѵ����Ԥ�����Ա�';['(R^2 =' num2str(R1) ' RMSE= ' num2str(error1) ' MSE= ' num2str(mse1) ' RPD= ' num2str(RPD1) ')' ]};
title(string)
%% Ԥ�⼯��ͼ
figure
plot(1:N,T_test,'r-*',1:N,T_sim2,'b-o','LineWidth',1.5)
legend('��ʵֵ','SSA-DBNԤ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ����')
string={'���Լ�Ԥ�����Ա�';['(R^2 =' num2str(R2) ' RMSE= ' num2str(error2)  ' MSE= ' num2str(mse2) ' RPD= ' num2str(RPD2) ')']};
title(string)

%% ���Լ����ͼ
figure  
ERROR3=T_test-T_sim2;
plot(T_test-T_sim2,'b-*','LineWidth',1.5)
xlabel('���Լ��������')
ylabel('Ԥ�����')
title('���Լ�Ԥ�����')
grid on;
legend('SSA-DBNԤ��������')
%% �����������ͼ
%% ѵ�������Ч��ͼ
figure
plot(T_train,T_sim1,'*r');
xlabel('��ʵֵ')
ylabel('Ԥ��ֵ')
string = {'ѵ����Ч��ͼ';['R^2_c=' num2str(R1)  '  RMSEC=' num2str(error1) ]};
title(string)
hold on ;h=lsline;
set(h,'LineWidth',1,'LineStyle','-','Color',[1 0 1])
%% Ԥ�⼯���Ч��ͼ
figure
plot(T_test,T_sim2,'ob');
xlabel('��ʵֵ')
ylabel('Ԥ��ֵ')
string1 = {'���Լ�Ч��ͼ';['R^2_p=' num2str(R2)  '  RMSEP=' num2str(error2) ]};
title(string1)
hold on ;h=lsline();
set(h,'LineWidth',1,'LineStyle','-','Color',[1 0 1])
%% ��ƽ��
R3=(R1+R2)./2;
error3=(error1+error2)./2;
%% ����������Ԥ�����ͼ
tsim=[T_sim1,T_sim2]';
S=[T_train,T_test]';
figure
plot(S,tsim,'ob');
xlabel('��ʵֵ')
ylabel('Ԥ��ֵ')
string1 = {'�����������Ԥ��ͼ';['R^2_p=' num2str(R3)  '  RMSEP=' num2str(error3) ]};
title(string1)
hold on ;h=lsline();
set(h,'LineWidth',1,'LineStyle','-','Color',[1 0 1])
%% ��ӡ������ָ��
disp(['-----------------------������--------------------------'])
disp(['ѵ�������۽��������ʾ��'])
disp(['ѵ����ƽ���������MAEΪ��',num2str(MAE1)])
disp(['ѵ�����������MSEΪ��       ',num2str(mse1)])
disp(['ѵ�������������RMSEΪ��  ',num2str(error1)])
disp(['ѵ��������ϵ��R^2Ϊ��  ',num2str(R1)])
disp(['ѵ����ʣ��Ԥ��в�RPDΪ��  ',num2str(RPD1)])
disp(['ѵ����ƽ�����԰ٷֱ����MAPEΪ��  ',num2str(MAPE1)])
disp(['���Լ����۽��������ʾ��'])
disp(['���Լ�ƽ���������MAEΪ��',num2str(MAE2)])
disp(['���Լ��������MSEΪ��       ',num2str(mse2)])
disp(['���Լ����������RMSEΪ��  ',num2str(error2)])
disp(['���Լ�����ϵ��R^2Ϊ��  ',num2str(R2)])
disp(['���Լ�ʣ��Ԥ��в�RPDΪ��  ',num2str(RPD2)])
disp(['���Լ�ƽ�����԰ٷֱ����MAPEΪ��  ',num2str(MAPE2)])
grid