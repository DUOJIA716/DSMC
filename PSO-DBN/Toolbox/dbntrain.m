function dbn = dbntrain(dbn, x, opts)

%%  训练深度置信网络
    n = numel(dbn.rbm);
    
%%  训练首层
    dbn.rbm{1} = rbmtrain(dbn.rbm{1}, x, opts);

%%  训练后续层
    for i = 2 : n
        x = rbmup(dbn.rbm{i - 1}, x);
        dbn.rbm{i} = rbmtrain(dbn.rbm{i}, x, opts);
    end

end