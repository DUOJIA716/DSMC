function labels = nnpredict(nn, x)

%%  ǰ�����
    nn = nnff(nn, x, zeros(size(x, 1), nn.size(end)));
    n = nn.n;

%%  �õ�������
    labels = nn.a{n};
    
end
