%% 颗粒碰撞处理
function [velocities, collisionDissipation, collisionCount] = handle_collisions(positions, velocities, params)
    collisionDissipation = 0;
    collisionCount = 0; % 初始化碰撞计数器
    grid = floor(positions / params.gridSize) + 1; % 网格划分
    uniqueCells = unique(grid, 'rows'); % 获取所有网格单元

    disp('Starting collision detection...');

    % 重力加速度
    g = params.gravity;

    for cellIdx = 1:size(uniqueCells, 1)
        % 获取当前网格内的颗粒
        cellParticles = ismember(grid, uniqueCells(cellIdx, :), 'rows');
        particlesInCell = find(cellParticles);

        % 检查填充比
        cellVolume = params.gridSize^3;
        particleVolume = 4/3 * pi * params.particleRadius^3;
        fillRatio = length(particlesInCell) * particleVolume / cellVolume;

        if fillRatio > params.fillRatioLimit
            warning(['Fill ratio exceeds limit in cell: ', mat2str(uniqueCells(cellIdx, :))]);
        end

        % 双重循环遍历颗粒对
        for i = 1:length(particlesInCell)-1
            for j = i+1:length(particlesInCell)
                p1 = particlesInCell(i);
                p2 = particlesInCell(j);

                relPos = positions(p1, :) - positions(p2, :);
                distance = norm(relPos);

                % 检查碰撞条件
                if distance < 2 * params.particleRadius
                    relVel = velocities(p1, :) - velocities(p2, :);
                    collisionNormal = relPos / max(distance, 1e-6); % 避免零向量
                    relSpeed = dot(relVel, collisionNormal);

                    if relSpeed < 0
                        % 保存碰撞前的速度
                        v1 = velocities(p1, :);
                        v2 = velocities(p2, :);

                        % 计算冲量
                        impulse = -params.restitutionCoeff * relSpeed * collisionNormal;

                        % 更新速度
                        velocities(p1, :) = v1 + 0.5 * impulse;
                        velocities(p2, :) = v2 - 0.5 * impulse;

                        % 碰撞后动能
                        postKE = 0.5 * params.particleMass * ...
                                 (norm(velocities(p1, :))^2 + norm(velocities(p2, :))^2);

                        % 碰撞后重力势能
                        postPotential = params.particleMass * g * (positions(p1, 2) + positions(p2, 2));

                        % 总能量
                        preTotalEnergy = 0.5 * params.particleMass * ...
                                         (norm(v1)^2 + norm(v2)^2) + ...
                                         params.particleMass * g * (positions(p1, 2) + positions(p2, 2));
                        postTotalEnergy = postKE + postPotential;

                        % 能量耗散
                        dissipation = preTotalEnergy - postTotalEnergy;
                        if dissipation > 0
                            collisionDissipation = collisionDissipation + dissipation;
                        end

                        % 增加碰撞计数器
                        collisionCount = collisionCount + 1;

                        % 调试输出
                        disp(['Collision detected between particles ', num2str(p1), ' and ', num2str(p2)]);
                    end
                end
            end
        end
    end

    % 打印总碰撞数
    disp(['Number of collisions detected: ', num2str(collisionCount)]);
    disp(['Collision dissipation: ', num2str(collisionDissipation)]);
end