function nn = nnsetup(architecture)

% NNSETUP creates a Feedforward Backpropagate Neural Network
% nn = nnsetup(architecture) returns an neural network structure with n=numel(architecture)
% layers, architecture being a n x 1 vector of layer sizes e.g. [784 100 10]

%%  �����������ĳ�ʼ��
    nn.size   = architecture;
    nn.n      = numel(nn.size);
    
    nn.activation_function              = 'tanh_opt';   %  �������'sigm' (sigmoid) or 'tanh_opt' (tanh).
    nn.learningRate                     = 2;            %  ѧϰ��
    nn.momentum                         = 0.5;          %  ��������
    nn.scaling_learningRate             = 1;            %  ѧϰ�ʵı�������
    nn.weightPenaltyL2                  = 0;            %  L2���򻯲���
    nn.nonSparsityPenalty               = 0;            %  ��ϡ��ͷ���
    nn.sparsityTarget                   = 0.05;         %  ϡ��Ŀ��
    nn.inputZeroMaskedFraction          = 0;            %  �Ƿ�ʹ��ȥ���Ա�����
    nn.dropoutFraction                  = 0;            %  dropout ��������
    nn.testing                          = 0;            %  �ڲ�������nntest ��������Ϊ 1
    nn.output                           = 'sigm';       %  ����������'sigm' (=logistic), 'softmax' and 'linear'

    for i = 2 : nn.n   
        % Ȩ�غ�Ȩ�ض���
        nn.W {i - 1} = (rand(nn.size(i), nn.size(i - 1) + 1) - 0.5) * 2 * 4 * sqrt(6 / (nn.size(i) + nn.size(i - 1)));
        nn.vW{i - 1} = zeros(size(nn.W{i - 1}));
        
        % ƽ���������ϡ���ԣ�
        nn.p{i} = zeros(1, nn.size(i));   
    end
end
