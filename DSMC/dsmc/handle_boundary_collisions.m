%% 边界碰撞处理
function [positions, velocities, boundaryDissipation] = handle_boundary_collisions(positions, velocities, params, structureVelocity)
    % 处理颗粒与外壳的边界碰撞
    % params: 模拟参数
    % structureVelocity: 外壳在z方向的速度

    boundaryDissipation = 0;

    % 场景1：颗粒与外壳底部和顶部之间的碰撞
    bottomCollision = positions(:, 3) < params.particleRadius;
    topCollision = positions(:, 3) > (params.boxSize(3) - params.particleRadius);

    if any(bottomCollision)
        % 底部碰撞速度更新（反弹 + 恢复系数）
        preVelocities = velocities(bottomCollision, 3); % 碰撞前速度
        velocities(bottomCollision, 3) = -params.restitutionCoeff * ...
                                         (velocities(bottomCollision, 3) - structureVelocity) + structureVelocity;

        % 能量耗散计算
        preKE = 0.5 * params.particleMass * preVelocities.^2;
        postKE = 0.5 * params.particleMass * velocities(bottomCollision, 3).^2;
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置（防止颗粒穿透边界）
        positions(bottomCollision, 3) = params.particleRadius;
    end

    if any(topCollision)
        % 顶部碰撞速度更新
        preVelocities = velocities(topCollision, 3);
        velocities(topCollision, 3) = -params.restitutionCoeff * velocities(topCollision, 3);

        % 能量耗散计算
        preKE = 0.5 * params.particleMass * preVelocities.^2;
        postKE = 0.5 * params.particleMass * velocities(topCollision, 3).^2;
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置
        positions(topCollision, 3) = params.boxSize(3) - params.particleRadius;
    end

    % 场景2：颗粒与外壳左右壁面和前后壁面的碰撞
    % 左右壁面碰撞
    leftCollision = positions(:, 1) < params.particleRadius;
    rightCollision = positions(:, 1) > (params.boxSize(1) - params.particleRadius);

    if any(leftCollision)
        % 左壁碰撞速度更新
        preVelocities = velocities(leftCollision, 1);
        velocities(leftCollision, 1) = -params.restitutionCoeff * velocities(leftCollision, 1);

        % 能量耗散计算
        preKE = 0.5 * params.particleMass * preVelocities.^2;
        postKE = 0.5 * params.particleMass * velocities(leftCollision, 1).^2;
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置
        positions(leftCollision, 1) = params.particleRadius;
    end

    if any(rightCollision)
        % 右壁碰撞速度更新
        preVelocities = velocities(rightCollision, 1);
        velocities(rightCollision, 1) = -params.restitutionCoeff * velocities(rightCollision, 1);

        % 能量耗散计算
        preKE = 0.5 * params.particleMass * preVelocities.^2;
        postKE = 0.5 * params.particleMass * velocities(rightCollision, 1).^2;
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置
        positions(rightCollision, 1) = params.boxSize(1) - params.particleRadius;
    end

    % 前后壁面碰撞
    frontCollision = positions(:, 2) < params.particleRadius;
    backCollision = positions(:, 2) > (params.boxSize(2) - params.particleRadius);

    if any(frontCollision)
        % 前壁碰撞速度更新
        preVelocities = velocities(frontCollision, 2);
        velocities(frontCollision, 2) = -params.restitutionCoeff * velocities(frontCollision, 2);

        % 能量耗散计算
        preKE = 0.5 * params.particleMass * preVelocities.^2;
        postKE = 0.5 * params.particleMass * velocities(frontCollision, 2).^2;
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置
        positions(frontCollision, 2) = params.particleRadius;
    end

    if any(backCollision)
        % 后壁碰撞速度更新
        preVelocities = velocities(backCollision, 2);
        velocities(backCollision, 2) = -params.restitutionCoeff * velocities(backCollision, 2);

        % 能量耗散计算
        preKE = 0.5 * params.particleMass * preVelocities.^2;
        postKE = 0.5 * params.particleMass * velocities(backCollision, 2).^2;
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置
        positions(backCollision, 2) = params.boxSize(2) - params.particleRadius;
    end
end

