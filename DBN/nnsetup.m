function nn = nnsetup(architecture)

% NNSETUP creates a Feedforward Backpropagate Neural Network
% nn = nnsetup(architecture) returns an neural network structure with n=numel(architecture)
% layers, architecture being a n x 1 vector of layer sizes e.g. [784 100 10]

%%  深度置信网络的初始化
    nn.size   = architecture;
    nn.n      = numel(nn.size);
    
    nn.activation_function              = 'tanh_opt';   %  激活函数：'sigm' (sigmoid) or 'tanh_opt' (tanh).
    nn.learningRate                     = 2;            %  学习率
    nn.momentum                         = 0.5;          %  动量参数
    nn.scaling_learningRate             = 1;            %  学习率的比例因子
    nn.weightPenaltyL2                  = 0;            %  L2正则化参数
    nn.nonSparsityPenalty               = 0;            %  非稀疏惩罚项
    nn.sparsityTarget                   = 0.05;         %  稀疏目标
    nn.inputZeroMaskedFraction          = 0;            %  是否使用去噪自编码器
    nn.dropoutFraction                  = 0;            %  dropout 参数因子
    nn.testing                          = 0;            %  内部变量，nntest 将此设置为 1
    nn.output                           = 'sigm';       %  输出激活函数，'sigm' (=logistic), 'softmax' and 'linear'

    for i = 2 : nn.n   
        % 权重和权重动量
        nn.W {i - 1} = (rand(nn.size(i), nn.size(i - 1) + 1) - 0.5) * 2 * 4 * sqrt(6 / (nn.size(i) + nn.size(i - 1)));
        nn.vW{i - 1} = zeros(size(nn.W{i - 1}));
        
        % 平均激活（用于稀疏性）
        nn.p{i} = zeros(1, nn.size(i));   
    end
end
