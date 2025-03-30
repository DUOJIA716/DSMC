%% 动态调整时间步长
function timeStep = adjust_time_step(velocities, params)
    % 根据最大速度动态调整时间步长
    maxSpeed = max(sqrt(sum(velocities.^2, 2)));
    timeStep = min(params.timeStep, params.gridSize / maxSpeed * params.safetyFactor);
end