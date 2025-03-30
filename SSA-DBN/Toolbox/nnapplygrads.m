function nn = nnapplygrads(nn)

% NNAPPLYGRADS updates weights and biases with calculated gradients
% nn = nnapplygrads(nn) returns an neural network structure with updated
% weights and biases

%%  权重移植
    % 权重导数
    for i = 1 : (nn.n - 1)
        
        % 判断是否存在正则化系数
        if(nn.weightPenaltyL2 > 0)
            dW = nn.dW{i} + nn.weightPenaltyL2 * nn.W{i};
        else
            dW = nn.dW{i};
        end
        
        % 更新权重
        dW = nn.learningRate * dW;
        
        % 判断是否存在动量
        if(nn.momentum > 0)
            nn.vW{i} = nn.momentum * nn.vW{i} + dW;
            dW = nn.vW{i};
        end
        
        % 更新权重
        nn.W{i} = nn.W{i} - dW;
        
    end
end
