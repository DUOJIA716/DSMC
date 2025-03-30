function nn = nnbp(nn)

% NNBP performs backpropagation
% nn = nnbp(nn) returns an neural network structure with updated weights 

%%  误差反向传播    
    n = nn.n;
    sparsityError = 0;

    % 输出层误差计算
    switch nn.output
        case 'sigm'
            d{n} = - nn.e .* (nn.a{n} .* (1 - nn.a{n}));
        case {'softmax', 'linear'}
            d{n} = - nn.e;
    end

    % 中间隐藏层误差计算
    for i = (n - 1) : -1 : 2
        % 激活函数的导数
        switch nn.activation_function 
            case 'sigm'
                d_act = nn.a{i} .* (1 - nn.a{i});
            case 'tanh_opt'
                d_act = 1.7159 * 2 / 3 * (1 - 1 / (1.7159)^2 * nn.a{i} .^ 2);
        end
        
        % 存在稀疏性
        if(nn.nonSparsityPenalty > 0)
            pi = repmat(nn.p{i}, size(nn.a{i}, 1), 1);
            sparsityError = [zeros(size(nn.a{i}, 1), 1), nn.nonSparsityPenalty ...
                * (-nn.sparsityTarget ./ pi + (1 - nn.sparsityTarget) ./ (1 - pi))];
        end
        
        % 反向传播一阶导数，是否删除偏置项
        if i + 1 == n
            d{i} = (d{i + 1} * nn.W{i} + sparsityError) .* d_act; 
        else
            d{i} = (d{i + 1}(:,2:end) * nn.W{i} + sparsityError) .* d_act;
        end
        
        % 是否存在dropout项
        if(nn.dropoutFraction > 0)
            d{i} = d{i} .* [ones(size(d{i}, 1), 1), nn.dropOutMask{i}];
        end

    end

    % 权重更新
    for i = 1 : (n - 1)
        if i + 1 == n
            nn.dW{i} = (d{i + 1}' * nn.a{i}) / size(d{i + 1}, 1);
        else
            nn.dW{i} = (d{i + 1}(:, 2: end)' * nn.a{i}) / size(d{i + 1}, 1);      
        end
    end

end