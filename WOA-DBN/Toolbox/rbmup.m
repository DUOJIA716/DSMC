function x = rbmup(rbm, x)

%%  ���޲��������� ���ϼ���
    x = sigm(repmat(rbm.c', size(x, 1), 1) + x * rbm.W');
    
end
