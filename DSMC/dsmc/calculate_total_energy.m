%% 计算总能量

function totalEnergy = calculate_total_energy(positions, velocities, params)
    % 计算系统的总能量
    kineticEnergy = 0.5 * params.particleMass * sum(sum(velocities.^2)); % 动能
    potentialEnergy = sum(params.particleMass * params.gravity * positions(:, 3)); % 势能
    totalEnergy = kineticEnergy + potentialEnergy;
end
