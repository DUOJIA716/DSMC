%% PSO-DBN
%% ��ջ�������
warning off             % �رվ�����Ϣ
close all               % �ر�����ͼ��
clear                   % ��չ���������
clc                     % ���������

%% ���·��
addpath('Toolbox\')      % ��ӹ�����·��

%% ��������
P_train = xlsread('data','training set','B2:G181')';  % ����ѵ������������
T_train = xlsread('data','training set','H2:H181')';  % ����ѵ�����������
P_test = xlsread('data','test set','B2:G45')';        % ������Լ���������
T_test = xlsread('data','test set','H2:H45')';        % ������Լ��������

%% ���ݲ���ע��
M = size(P_train, 2);  % ѵ��������
N = size(P_test, 2);   % ����������
outdim = 1;            % ���ά��

%% ���ݹ�һ��
[p_train, ps_input] = mapminmax(P_train, 0, 1);    % �������ݹ�һ��
p_test = mapminmax('apply', P_test, ps_input);     % Ӧ�ù�һ�����������Լ���������
[t_train, ps_output] = mapminmax(T_train, 0, 1);   % ������ݹ�һ��
t_test = mapminmax('apply', T_test, ps_output);    % Ӧ�ù�һ�����������Լ��������

%% ����ת������Ӧģ��
p_train = p_train'; 
p_test = p_test';
t_train = t_train'; 
t_test = t_test';

%% ��ʼ��opts����
opts.numepochs = 500;           % ѵ������
opts.batchsize = 12;            % ÿ��ѵ����������
opts.momentum  = 0;             % ��������
opts.alpha     = 0.05;          % ѧϰ��

%% PSO�Ż��Ĳ�������
dim = 4;                          % �Ż�������ά��
lb  = [20, 20, 1000, 1.5];        % �Ż���������
ub  = [80, 80, 3000, 3.0];        % �Ż���������
numsum = length(ub) - 2;          % ���ز�ڵ���
fun = @(x)fun(x, numsum, p_train, t_train, opts);  % ����Ŀ�꺯��

pop = 6;                           % ����Ⱥ��ģ
Max_iteration = 6;                 % ����������
Vmax = 0.1;                        % �ٶ�����
Vmin = -0.1;                       % �ٶ�����

%% PSO�㷨�Ż�DBN����
tic;  % ��Ӽ�ʱ��ʼ
[Best_pos, Best_score, curve] = PSO(pop, Max_iteration, lb, ub, dim, fun, Vmax, Vmin); 
zbest = Best_pos;                 % ��ȡ���λ��

%% ��ӡ���Ų���
disp('------------------���Ų������------------------')
disp(['���ز�ڵ���: ', num2str(round(zbest(1:numsum)))])
disp(['ѧϰ��: ', num2str(zbest(end))])
disp(['��������: ', num2str(round(zbest(numsum + 1)))])
disp(['������Ӧ��ֵ: ', num2str(Best_score)])
disp('------------------------------------------------')

%% ��ȡ���Ų���
for i = 1:numsum
    dbn.sizes(i + 1) = round(zbest(i));  % ����ÿ��Ľڵ���
end
dbn.sizes(1) = [];                       % ɾ�������ڵ���
lr = zbest(end);                         % ��ȡѧϰ��
epochs = round(zbest(numsum + 1));       % ��ȡ��������

%% ѵ��DBNģ��
dbn = dbnsetup(dbn, p_train, opts);      % ����DBNģ��
dbn = dbntrain(dbn, p_train, opts);      % ѵ��DBNģ��

%% �������㲢����΢��
nn = dbnunfoldtonn(dbn, outdim);         % չ��DBN��������
nn.activation_function = 'sigm';         % ���ü����
nn.learningRate = lr;                    % ����ѧϰ��
nn.momentum = 0.5;                       % ���ö�������

[nn, loss] = nntrain(nn, p_train, t_train, opts);  % ����΢��ѵ��

%% Ԥ����
t_sim1 = nnpredict(nn, p_train);         % ѵ����Ԥ��
t_sim2 = nnpredict(nn, p_test);          % ���Լ�Ԥ��

%%  �������ʱ�䲢��ʾ
simulationTime = toc;  % ��ʱ����
disp(['��ֵ�����ܺ�ʱ: ', num2str(simulationTime), ' ��']);

%% ���ݷ���һ��
T_sim1 = mapminmax('reverse', t_sim1', ps_output);  % ѵ��������һ��
T_sim2 = mapminmax('reverse', t_sim2', ps_output);  % ���Լ�����һ��

%% ������ʧ��������
figure
plot(1:length(loss), loss, 'b-', 'LineWidth', 1)
xlim([1, length(loss)])
xlabel('��������')
ylabel('�����ʧ')
legend('��ʧ����')
title('��ʧ��������')
grid on

%% ������������
figure
plot(curve, 'linewidth', 1.5);
grid on
xlabel('��������')
ylabel('��Ӧ��ֵ')
title('PSO-DBN��������')

%% ���Լ��ع�ͼ��������
figure;
plotregression(T_test, T_sim2, '�ع�ͼ');  % ���ƻع�ͼ
figure;
ploterrhist(T_test - T_sim2, '���ֱ��ͼ');  % �������ֱ��ͼ

%% ��������ָ��
error1 = sqrt(sum((T_sim1 - T_train).^2) ./ M);  % ѵ�������������RMSE
error2 = sqrt(sum((T_test - T_sim2).^2) ./ N);  % ���Լ����������RMSE
R1 = 1 - norm(T_train - T_sim1)^2 / norm(T_train - mean(T_train))^2;  % ѵ��������ϵ��
R2 = 1 - norm(T_test - T_sim2)^2 / norm(T_test - mean(T_test))^2;     % ���Լ�����ϵ��
%������� MSE
mse1 = sum((T_sim1 - T_train).^2)./M;
mse2 = sum((T_sim2 - T_test).^2)./N;

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

%% ����Ԥ�����Ա�ͼ
figure
plot(1:M, T_train, 'r-*', 1:M, T_sim1, 'b-o', 'LineWidth', 1.5)
legend('��ʵֵ', 'PSO-DBNԤ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ����')
title(['ѵ����Ԥ�����Ա� (R^2 =', num2str(R1), ', RMSE=', num2str(error1), ')'])

figure
plot(1:N, T_test, 'r-*', 1:N, T_sim2, 'b-o', 'LineWidth', 1.5)
legend('��ʵֵ', 'PSO-DBNԤ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ����')
title(['���Լ�Ԥ�����Ա� (R^2 =', num2str(R2), ', RMSE=', num2str(error2), ')'])

%% ���Լ����ͼ
figure  
ERROR3 = T_test - T_sim2;  % �������
plot(ERROR3, 'b-*', 'LineWidth', 1.5)
xlabel('���Լ��������')
ylabel('Ԥ�����')
title('���Լ�Ԥ�����')
grid on;
legend('PSO-DBNԤ��������')

%% ��ӡ����ָ��
disp('-------------------���Լ�����ָ��-------------------')
disp(['���������RMSE: ', num2str(error2)])
disp(['����ϵ��R^2: ', num2str(R2)])


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

