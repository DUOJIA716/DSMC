function [nn, L]  = nntrain(nn, train_x, train_y, opts)

% NNTRAIN trains a neural net
% [nn, L] = nnff(nn, x, y, opts) trains the neural network nn with input x and
% output y for opts.numepochs epochs, with minibatches of size
% opts.batchsize. Returns a neural network nn with updated activations,
% errors, weights and biases, (nn.a, nn.e, nn.W, nn.b) and L, the sum
% squared error for each training minibatch.

%%  �����������ͺ͸���
assert(isfloat(train_x), 'train_x must be a float');
assert(nargin == 4, 'number ofinput arguments must be 4')

%%  ÿ��ѵ��������������Ϊ�������������ܱ�ȫ����������M����
M = size(train_x, 1);
batchsize = opts.batchsize;
numepochs = opts.numepochs;
numbatches = M / batchsize;
assert(rem(numbatches, 1) == 0, 'numbatches must be a integer');

%%  ��ʧ��������
L = zeros(numepochs, 1);

%%  ģ�ͷ���΢��
for i = 1 : numepochs
    % Ԥ���� 
    batch_loss = 0;
    kk = randperm(M);
    for l = 1 : numbatches
        batch_x = train_x(kk((l - 1) * batchsize + 1 : l * batchsize), :);
        
        % Add noise to input (for use in denoising autoencoder)
        % ��������������У�ȥ���Ա�����
        if(nn.inputZeroMaskedFraction ~= 0)
            batch_x = batch_x .* (rand(size(batch_x)) > nn.inputZeroMaskedFraction);
        end
        
        batch_y = train_y(kk((l - 1) * batchsize + 1 : l * batchsize), :);
        
        % ǰ����㣬����΢����Ȩ����ֲ
        nn = nnff(nn, batch_x, batch_y);
        nn = nnbp(nn);
        nn = nnapplygrads(nn);
        
        % �õ�ƽ����ʧ
        batch_loss = batch_loss + nn.L / numbatches;
        
    end
    
    % ����ѧϰ��
    nn.learningRate = nn.learningRate * nn.scaling_learningRate;

    % ��¼��ʧ��������
    L(i) = batch_loss;
end

end