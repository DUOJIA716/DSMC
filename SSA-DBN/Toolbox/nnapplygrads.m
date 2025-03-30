function nn = nnapplygrads(nn)

% NNAPPLYGRADS updates weights and biases with calculated gradients
% nn = nnapplygrads(nn) returns an neural network structure with updated
% weights and biases

%%  Ȩ����ֲ
    % Ȩ�ص���
    for i = 1 : (nn.n - 1)
        
        % �ж��Ƿ��������ϵ��
        if(nn.weightPenaltyL2 > 0)
            dW = nn.dW{i} + nn.weightPenaltyL2 * nn.W{i};
        else
            dW = nn.dW{i};
        end
        
        % ����Ȩ��
        dW = nn.learningRate * dW;
        
        % �ж��Ƿ���ڶ���
        if(nn.momentum > 0)
            nn.vW{i} = nn.momentum * nn.vW{i} + dW;
            dW = nn.vW{i};
        end
        
        % ����Ȩ��
        nn.W{i} = nn.W{i} - dW;
        
    end
end
