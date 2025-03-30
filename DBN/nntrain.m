function [nn, L]  = nntrain(nn, train_x, train_y, opts)

% NNTRAIN trains a neural net
% [nn, L] = nnff(nn, x, y, opts) trains the neural network nn with input x and
% output y for opts.numepochs epochs, with minibatches of size
% opts.batchsize. Returns a neural network nn with updated activations,
% errors, weights and biases, (nn.a, nn.e, nn.W, nn.b) and L, the sum
% squared error for each training minibatch.

%%  限制输入类型和个数
assert(isfloat(train_x), 'train_x must be a float');
assert(nargin == 4, 'number ofinput arguments must be 4')

%%  每次训练样本个数必须为整数个数，且能被全部样本个数M整除
M = size(train_x, 1);
batchsize = opts.batchsize;
numepochs = opts.numepochs;
numbatches = M / batchsize;
assert(rem(numbatches, 1) == 0, 'numbatches must be a integer');

%%  损失函数曲线
L = zeros(numepochs, 1);

%%  模型反向微调
for i = 1 : numepochs
    % 预定义 
    batch_loss = 0;
    kk = randperm(M);
    for l = 1 : numbatches
        batch_x = train_x(kk((l - 1) * batchsize + 1 : l * batchsize), :);
        
        % Add noise to input (for use in denoising autoencoder)
        % 添加噪声到输入中，去噪自编码器
        if(nn.inputZeroMaskedFraction ~= 0)
            batch_x = batch_x .* (rand(size(batch_x)) > nn.inputZeroMaskedFraction);
        end
        
        batch_y = train_y(kk((l - 1) * batchsize + 1 : l * batchsize), :);
        
        % 前向计算，反向微调，权重移植
        nn = nnff(nn, batch_x, batch_y);
        nn = nnbp(nn);
        nn = nnapplygrads(nn);
        
        % 得到平均损失
        batch_loss = batch_loss + nn.L / numbatches;
        
    end
    
    % 更新学习率
    nn.learningRate = nn.learningRate * nn.scaling_learningRate;

    % 记录损失函数曲线
    L(i) = batch_loss;
end

end