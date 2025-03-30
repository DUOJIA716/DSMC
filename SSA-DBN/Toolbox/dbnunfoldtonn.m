function nn = dbnunfoldtonn(dbn, outputsize)

%%  结构初始化，添加输出层
    if(exist('outputsize', 'var'))
        size = [dbn.sizes, outputsize];
    else
        size = [dbn.sizes];
    end

%%  初始化模型
    nn = nnsetup(size);

%%  权重抑制
    for i = 1 : numel(dbn.rbm)
        nn.W{i} = [dbn.rbm{i}.c, dbn.rbm{i}.W];
    end
    
end