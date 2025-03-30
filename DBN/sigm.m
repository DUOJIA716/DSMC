function X = sigm(P)

%% 激活函数
X = 1 ./ (1 + exp(-P));

end