function X = sigmrnd(P)

%%  ¼¤»î¸ÅÂÊ
    X = double(1 ./ (1 + exp(-P)) > rand(size(P)));
    
end