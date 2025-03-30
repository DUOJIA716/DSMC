function fitness = fun(pop, numsum, p_train, t_train, opts)

%%  ��ȡ����
for i = 1 : numsum
    dbn.sizes(i + 1) = round(pop(i)); % ���ز�ڵ�
end

%%  ��ȡ���Ų���
lr = pop(end);                        % ѧϰ��
epochs = round(pop(numsum + 1));      % ��������
dbn.sizes(1) = [];

%%  ���ݵĲ���
num_size = length(t_train);

%%  ������֤����
indices = crossvalind('Kfold', num_size, 3);

for i = 1 : 3
    
    % ��ȡ��i�����ݵ������߼�ֵ
    valid_data = (indices == i);
    
    % ȡ������ȡ��i��ѵ�����ݵ������߼�ֵ
    train_data = ~valid_data;
    
    % 1�ݲ��ԣ�2��ѵ��
    pv_train = p_train(train_data, :);
    tv_train = t_train(train_data, :);
    
    pv_valid = p_train(valid_data, :);
    tv_valid = t_train(valid_data, :);
    
    
    % Ԥѵ���Ĳ�������
    Ms = size(pv_train , 1);                % ѵ����������
    opts.numepochs = 300;                   % ѵ������
    opts.batchsize = 12;                    % ÿ��ѵ���������� �����㣺��Ms / batchsize = ������
    opts.momentum  = 0;                     % ѧϰ�ʵĶ���
    opts.alpha     = 0.05;                  % ѧϰ��
    
    % ѵ��ģ��
    dbn = dbnsetup(dbn, pv_train, opts);     % ����ģ��
    dbn = dbntrain(dbn, pv_train, opts);     % ѵ��ģ��

    % ѵ��Ȩ����ֲ����������
    outdim = size(tv_train, 2);
    nn = dbnunfoldtonn(dbn, outdim);

    % �����������
    opts.numepochs          = epochs;       % ����΢������
    opts.batchsize          = 12;           % ÿ�η���΢�������� �����㣺��Ms / batchsize = ������
    
    nn.activation_function  = 'sigm';       % �����
    nn.learningRate         = lr;           % ѧϰ��
    nn.momentum             = 0.5;          % ��������
    nn.scaling_learningRate = 1;            % ѧϰ�ʵı�������
    
    [net, ~] = nntrain(nn, pv_train, tv_train, opts);  % ѵ��

    % Ԥ�� 
    t_sim = nnpredict(net, pv_valid); 

    % �õ���Ӧ��ֵ
    error(i) = mean(sqrt(sum((t_sim - tv_valid) .^ 2, 2) ./ size(pv_valid, 1)));

end

%%  ��ȡ��Ӧ��
fitness = mean(error);