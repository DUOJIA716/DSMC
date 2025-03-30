function labels = nnpredict(nn, x)

%%  前向计算
    nn = nnff(nn, x, zeros(size(x, 1), nn.size(end)));
    n = nn.n;

%%  得到输出结果
    labels = nn.a{n};
    
end
