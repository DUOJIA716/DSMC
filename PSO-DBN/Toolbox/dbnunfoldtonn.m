function nn = dbnunfoldtonn(dbn, outputsize)

%%  �ṹ��ʼ������������
    if(exist('outputsize', 'var'))
        size = [dbn.sizes, outputsize];
    else
        size = [dbn.sizes];
    end

%%  ��ʼ��ģ��
    nn = nnsetup(size);

%%  Ȩ������
    for i = 1 : numel(dbn.rbm)
        nn.W{i} = [dbn.rbm{i}.c, dbn.rbm{i}.W];
    end
    
end