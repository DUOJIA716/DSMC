function x = rbmup(rbm, x)

%%  受限玻尔兹曼机 向上计算
    x = sigm(repmat(rbm.c', size(x, 1), 1) + x * rbm.W');
    
end
