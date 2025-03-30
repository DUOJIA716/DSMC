%% 分析和可视化结果

function analyze_results(totalEnergyHistory, positions, params)
    % 能量分析
    figure;
    plot(totalEnergyHistory);
    xlabel('时间步');
    ylabel('总能量 (J)');
    title('总能量随时间的变化');

    % 最终颗粒分布
    figure;
    scatter3(positions(:, 1), positions(:, 2), positions(:, 3), 20, 'filled');
    xlim([0, params.boxSize(1)]);
    ylim([0, params.boxSize(2)]);
    zlim([0, params.boxSize(3)]);
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('最终颗粒分布');
    grid on;


end