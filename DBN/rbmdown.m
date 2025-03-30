function x = rbmdown(rbm, x)

%%  受限玻尔兹曼机 向下计算
    x = sigm(repmat(rbm.b', size(x, 1), 1) + x * rbm.W);
end
