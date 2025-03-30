function fitness = fun(pop, numsum, p_train, t_train, opts)

%%  提取参数
for i = 1 : numsum
    dbn.sizes(i + 1) = round(pop(i)); % 隐藏层节点
end

%%  提取最优参数
lr = pop(end);                        % 学习率
epochs = round(pop(numsum + 1));      % 迭代次数
dbn.sizes(1) = [];

%%  数据的参数
num_size = length(t_train);

%%  交叉验证程序
indices = crossvalind('Kfold', num_size, 3);

for i = 1 : 3
    
    % 获取第i份数据的索引逻辑值
    valid_data = (indices == i);
    
    % 取反，获取第i份训练数据的索引逻辑值
    train_data = ~valid_data;
    
    % 1份测试，2份训练
    pv_train = p_train(train_data, :);
    tv_train = t_train(train_data, :);
    
    pv_valid = p_train(valid_data, :);
    tv_valid = t_train(valid_data, :);
    
    
    % 预训练的参数设置
    Ms = size(pv_train , 1);                % 训练集样本数
    opts.numepochs = 300;                   % 训练次数
    opts.batchsize = 12;                    % 每次训练样本个数 需满足：（Ms / batchsize = 整数）
    opts.momentum  = 0;                     % 学习率的动量
    opts.alpha     = 0.05;                  % 学习率
    
    % 训练模型
    dbn = dbnsetup(dbn, pv_train, opts);     % 建立模型
    dbn = dbntrain(dbn, pv_train, opts);     % 训练模型

    % 训练权重移植，添加输出层
    outdim = size(tv_train, 2);
    nn = dbnunfoldtonn(dbn, outdim);

    % 反向调整网络
    opts.numepochs          = epochs;       % 反向微调次数
    opts.batchsize          = 12;           % 每次反向微调样本数 需满足：（Ms / batchsize = 整数）
    
    nn.activation_function  = 'sigm';       % 激活函数
    nn.learningRate         = lr;           % 学习率
    nn.momentum             = 0.5;          % 动量参数
    nn.scaling_learningRate = 1;            % 学习率的比例因子
    
    [net, ~] = nntrain(nn, pv_train, tv_train, opts);  % 训练

    % 预测 
    t_sim = nnpredict(net, pv_valid); 

    % 得到适应度值
    error(i) = mean(sqrt(sum((t_sim - tv_valid) .^ 2, 2) ./ size(pv_valid, 1)));

end

%%  获取适应度
fitness = mean(error);