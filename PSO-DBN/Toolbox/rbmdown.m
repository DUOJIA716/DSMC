function x = rbmdown(rbm, x)

%%  ���޲��������� ���¼���
    x = sigm(repmat(rbm.b', size(x, 1), 1) + x * rbm.W);
end
