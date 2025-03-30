%% 初始化颗粒
function [positions, velocities] = initialize_particles(params)
    % 初始化颗粒的位置和速度

    positions = [rand(params.numParticles, 2) .* (params.boxSize(1:2) - 2 * params.particleRadius), ...
             params.particleRadius + rand(params.numParticles, 1) * (params.boxSize(3) - 2 * params.particleRadius)];


    velocities = zeros(params.numParticles, 3); % 初始速度设置为零
    % velocities = (rand(params.numParticles, 3) - 0.5) * 0.2; % 实验颗粒速度
end
