function dbn = dbntrain(dbn, x, opts)

%%  ѵ�������������
    n = numel(dbn.rbm);
    
%%  ѵ���ײ�
    dbn.rbm{1} = rbmtrain(dbn.rbm{1}, x, opts);

%%  ѵ��������
    for i = 2 : n
        x = rbmup(dbn.rbm{i - 1}, x);
        dbn.rbm{i} = rbmtrain(dbn.rbm{i}, x, opts);
    end

end