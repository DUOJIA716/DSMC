%% 更新颗粒位置
function positions = update_positions(positions, velocities, params)
    % 根据速度更新颗粒位置
    positions = positions + velocities * params.timeStep;

end