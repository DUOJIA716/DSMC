%% 应用重力
function velocities = apply_gravity(velocities, params)
    % 为颗粒速度施加重力加速度
    velocities(:, 3) = velocities(:, 3) - params.gravity * params.timeStep;
end