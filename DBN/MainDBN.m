%% DBN
%% ��ջ�������
warning off             % �رձ�����Ϣ
close all               % �رտ�����ͼ��
clear                   % ��ձ���
clc                     % ���������

%% ��������
P_train = xlsread('data', 'training set', 'B2:G181')';  % ����ѵ������������
T_train = xlsread('data', 'training set', 'H2:H181')';  % ����ѵ�����������
P_test  = xlsread('data', 'test set', 'B2:G45')';       % ������Լ���������
T_test  = xlsread('data', 'test set', 'H2:H45')';       % ������Լ��������

M = size(P_train, 2);   % ѵ������������
N = size(P_test, 2);    % ���Լ���������

%% ���ݹ�һ��
[p_train, ps_input] = mapminmax(P_train, 0, 1);
p_test = mapminmax('apply', P_test, ps_input);

[t_train, ps_output] = mapminmax(T_train, 0, 1);
t_test = mapminmax('apply', T_test, ps_output);

%% ת������Ӧģ��
p_train = p_train'; p_test = p_test';
t_train = t_train'; t_test = t_test';

%% ģ��Ԥѵ��
dbn.sizes      = [50, 50];           % ���ز�ڵ���
opts.numepochs = 500;                % ѵ������
opts.batchsize = 12;                 % ÿ��ѵ�����������������㣺��M / batchsize = ������
opts.momentum  = 0;                  % ��������
opts.alpha     = 0.05;                % ѧϰ��

dbn = dbnsetup(dbn, p_train, opts);  % ����ģ��
dbn = dbntrain(dbn, p_train, opts);  % ѵ��ģ��

%% ѵ��Ȩ����ֲ����������
tic;  % ��Ӽ�ʱ��ʼ
nn = dbnunfoldtonn(dbn, size(T_train, 1));

%% �����������
opts.numepochs = 2000;               % ����΢������
opts.batchsize = 12;                 % ÿ�η���΢���������������㣺��M / batchsize = ������

nn.activation_function = 'sigm';     % �����
nn.learningRate = 0.4;                 % ѧϰ��
nn.momentum = 0.5;                   % ��������
nn.scaling_learningRate = 1;         % ѧϰ�ʵı�������

[nn, loss] = nntrain(nn, p_train, t_train, opts); % ����΢��ѵ��

%% ģ��Ԥ��
t_sim1 = nnpredict(nn, p_train);
t_sim2 = nnpredict(nn, p_test);

%%  �������ʱ�䲢��ʾ
simulationTime = toc;  % ��ʱ����
disp(['��ֵ�����ܺ�ʱ: ', num2str(simulationTime), ' ��']);

%% ���ݷ���һ��
T_sim1 = mapminmax('reverse', t_sim1, ps_output);
T_sim2 = mapminmax('reverse', t_sim2, ps_output);

%% ���Լ�����ع�ͼ��������
figure;
plotregression(T_test, T_sim2, '�ع�ͼ');
figure;
ploterrhist(T_test - T_sim2, '���ֱ��ͼ');

%% ���������
error1 = sqrt(sum((T_sim1' - T_train).^2) ./ M);
error2 = sqrt(sum((T_sim2' - T_test ).^2) ./ N);

%% ������ʧ��������
figure
plot(1:length(loss), loss, 'b-', 'LineWidth', 1)
xlim([1, length(loss)])
xlabel('��������')
ylabel('�����ʧ')
legend('��ʧ����')
title('��ʧ����')
grid

%% ��ͼ
figure
plot(1:M, T_train, '-s', 1:M, T_sim1, '-o', 'LineWidth', 1)
legend('��ʵֵ', 'Ԥ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ����')
string = {'ѵ����Ԥ�����Ա�'; ['RMSE=' num2str(error1)]};
title(string)
xlim([1, M])
grid

figure
plot(1:N, T_test, '-s', 1:N, T_sim2, '-o', 'LineWidth', 1)
legend('��ʵֵ', 'Ԥ��ֵ')
xlabel('Ԥ������')
ylabel('Ԥ����')
string = {'���Լ�Ԥ�����Ա�'; ['RMSE=' num2str(error2)]};
title(string)
xlim([1, N])
grid

%% ��������ָ��
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
% RPD��ʣ��Ԥ��в
rpd1 = std(T_train) / std(T_train - T_sim1');
rpd2 = std(T_test) / std(T_test - T_sim2');

%% �����������ͼ
% ѵ�������Ч��ͼ
figure
plot(T_train, T_sim1, '*r');
xlabel('��ʵֵ');
ylabel('Ԥ��ֵ');
string = {'ѵ����Ч��ͼ'; ['R^2_c=' num2str(R1) '  RMSEC=' num2str(error1)]};
title(string);
hold on;
h = lsline; % ������������
set(h, 'LineWidth', 1, 'LineStyle', '-', 'Color', [1 0 1]); % �����������ʽ

% ���Լ����Ч��ͼ
figure
plot(T_test, T_sim2, 'ob');
xlabel('��ʵֵ');
ylabel('Ԥ��ֵ');
string1 = {'���Լ�Ч��ͼ'; ['R^2_p=' num2str(R2) '  RMSEP=' num2str(error2)]};
title(string1);
hold on;
h = lsline(); % ������������
set(h, 'LineWidth', 1, 'LineStyle', '-', 'Color', [1 0 1]); % �����������ʽ


%% ��ӡ����ָ��
disp('-------------------����ָ��-------------------')
disp(['ѵ���� R^2: ', num2str(R1)])
disp(['���Լ� R^2: ', num2str(R2)])
disp(['ѵ���� MSE: ', num2str(mse1)])
disp(['���Լ� MSE: ', num2str(mse2)])
disp(['ѵ���� MAE: ', num2str(mae1)])
disp(['���Լ� MAE: ', num2str(mae2)])
disp(['ѵ���� MAPE: ', num2str(mape1), '%'])
disp(['���Լ� MAPE: ', num2str(mape2), '%'])
disp(['ѵ���� RMSE: ', num2str(rmse1)])
disp(['���Լ� RMSE: ', num2str(rmse2)])
disp(['ѵ���� RPD: ', num2str(rpd1)])
disp(['���Լ� RPD: ', num2str(rpd2)])
disp('---------------------------------------------')


