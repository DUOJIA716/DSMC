function X = sigmrnd(P)

%%  激活概率
    X = double(1 ./ (1 + exp(-P)) > rand(size(P)));
    
end