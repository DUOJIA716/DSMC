%% 更新结构动力学
function [newDisplacement, newVelocity] = update_structure(externalForce, displacement, velocity, params)
    % 更新结构的位移和速度（基于牛顿第二定律）
    acceleration = (externalForce - params.damping * velocity - params.stiffness * displacement) ...
        / params.massStructure;

    % 数值积分更新速度和位移
    newVelocity = velocity + acceleration * params.timeStep;
    newDisplacement = displacement + newVelocity * params.timeStep;
end