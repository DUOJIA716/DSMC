function nn = nnff(nn, x, y)

% NNFF performs a feedforward pass
% nn = nnff(nn, x, y) returns an neural network structure with updated
% layer activations, error and loss (nn.a, nn.e and nn.L)

%%  前向计算
    n = nn.n;
    m = size(x, 1);

    %%  网络第一层输出
    x = [ones(m, 1), x];
    nn.a{1} = x;

    %%  网络中间层的计算
    for i = 2 : n - 1

        % 选择激活函数，计算单元的输出
        switch nn.activation_function 
            case 'sigm'
                nn.a{i} = sigm(nn.a{i - 1} * nn.W{i - 1}');
            case 'tanh_opt'
                nn.a{i} = tanh_opt(nn.a{i - 1} * nn.W{i - 1}');
        end
        
        % 是否存在dropout层
        if(nn.dropoutFraction > 0)
            if(nn.testing)
                nn.a{i} = nn.a{i} .* (1 - nn.dropoutFraction);
            else
                nn.dropOutMask{i} = (rand(size(nn.a{i})) > nn.dropoutFraction);
                nn.a{i} = nn.a{i} .* nn.dropOutMask{i};
            end
        end
        
        % 计算用于稀疏性的运行指数激活
        if(nn.nonSparsityPenalty > 0)
            nn.p{i} = 0.99 * nn.p{i} + 0.01 * mean(nn.a{i}, 1);
        end
        
        % 添加偏置项
        nn.a{i} = [ones(m, 1), nn.a{i}];
    end

    % 选择输出层激活函数，得到最后一层输出
    switch nn.output 
        case 'sigm'
            nn.a{n} = sigm(nn.a{n - 1} * nn.W{n - 1}');
        case 'linear'
            nn.a{n} = nn.a{n - 1} * nn.W{n - 1}';
        case 'softmax'
            nn.a{n} = nn.a{n - 1} * nn.W{n - 1}';
            nn.a{n} = exp(bsxfun(@minus, nn.a{n}, max(nn.a{n},[],2)));
            nn.a{n} = bsxfun(@rdivide, nn.a{n}, sum(nn.a{n}, 2)); 
    end

    % 得到误差值
    nn.e = y - nn.a{n};
    
    % 选择激活函数得到损失函数值
    switch nn.output
        case {'sigm', 'linear'}
            nn.L = 1 / 2 * sum(sum(nn.e .^ 2)) / m; 
        case 'softmax'
            nn.L = -sum(sum(y .* log(nn.a{n}))) / m;
    end

end