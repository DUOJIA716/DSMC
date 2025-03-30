% 主程序：基于DSMC方法的颗粒阻尼模拟
clear; clc;

% 参数定义

% 主结构参数
params.structure_mass = 0.0376;                             % 主结构质量 (kg)
params.structure_stiffness = 470;                           % 刚度系数 K (N/m)
params.structure_damping = 0.1;                             % 阻尼系数 C (N·s/m)
params.structure_diameter = 25.4e-3;                        % 容器直径 (m)
params.structure_height = 7.04e-3;                          % 容器高度 (m)
params.structure_young_modulus = 7e10;                      % 弹性模量 (Pa)
params.structure_nu = 0.34;                                 % 泊松比
params.structure_coefficient_restitution = 0.5;             % 恢复系数
params.structure_coefficient_static_friction = 0.154;       % 静摩擦系数
params.structure_coefficient_rolling_friction = 0.01;       % 滚动摩擦系数
params.Amplitude = 0.0005; % 外部荷载振幅，单位：m
params.forceAmplitude = params.Amplitude * params.structure_stiffness; % 外部荷载，单位：N
params.forceFrequency = 15; % 外部荷载频率，单位：Hz

% 颗粒参数
params.particle_density = 1.13e4;                           % 颗粒密度 (kg/m^3)
params.particle_diameter = 0.88e-3;                         % 颗粒直径 (m)
params.particle_radius = params.particle_diameter / 2;
params.particle_mass = (4/3) * pi * (params.particle_radius^3) * params.particle_density;
params.particle_young_modulus = 1.6e10;                     % 弹性模量 (Pa)
params.particle_nu = 0.35;                                  % 泊松比
params.particle_coefficient_restitution = 0.7;              % 恢复系数
params.particle_coefficient_static_friction = 0.2;          % 静摩擦系数
params.particle_coefficient_rolling_friction = 0.01;        % 滚动摩擦系数
params.particle_number = 300;                              % 颗粒数量

% 通用参数
params.gravity = 9.81;                                  % 重力加速度 (m/s^2)
params.timeStep = 1e-5;                          % 初始时间步长 (s)
params.numSteps = 500000; % 总时间步数
params.total_time = params.timeStep * params.numSteps;      % 总模拟时间 (s)
params.safetyFactor = 0.8; % 时间步长调整的安全系数
params.fillRatioLimit = 0.63; % 最大填充比
params.gridSize = params.particle_diameter * 5; % 网格单元边长，单位：米

% 初始化颗粒
[positions, velocities] = initialize_particles(params);
totalparticleEnergyHistory = zeros(params.numSteps, 1); % 保存颗粒能量的历史数据

% 初始化主结构运动
structureDisplacement = 0; % 主结构初始位移
structureVelocity = 0; % 主结构初始速度
structureDisplacementHistory = zeros(1, params.numSteps); % 记录主系统位移
structureVelocityHistory = zeros(1, params.numSteps); % 记录主系统速度

% 初始化反作用力为零
particlesreactionForce = [0, 0, 0]; % 第一次迭代时使用零反作用力

% 初始化耗能记录
totalBoundaryDissipation = zeros(1, params.numSteps); % 累积边界碰撞耗能
totalCollisionDissipation = zeros(1, params.numSteps); % 累积颗粒间碰撞耗能

% 初始化全局颗粒碰撞次数
totalCollisionCount = 0;

% 模拟循环
for step = 1:params.numSteps
    % 当前时间
    currentTime = step * params.timeStep;

     % 计算外部荷载
    externalForce = params.forceAmplitude * sin(2 * pi * params.forceFrequency * currentTime);

     % 更新结构动力学
    [structureDisplacement, structureVelocity] = update_structure(externalForce, structureDisplacement, structureVelocity, particlesreactionForce, params);

    % 记录主系统位移和速度
    structureDisplacementHistory(step) = structureDisplacement;
    structureVelocityHistory(step) = structureVelocity;

    % 自由运动阶段
    positions = update_positions(positions, velocities, params);
    velocities = apply_gravity(velocities, params);

    % 动态调整时间步长
    params.timeStep = adjust_time_step(velocities, params);

    % 边界碰撞处理阶段
    [positions, velocities, boundaryDissipation, particlesreactionForce] = handle_boundary_collisions(positions, velocities, params, structureVelocity);
    totalBoundaryDissipation(step) = totalBoundaryDissipation(max(step - 1, 1)) + boundaryDissipation;

    % 颗粒碰撞处理阶段
    [velocities, collisionDissipation, collisionCount] = handle_collisions(positions, velocities, params);
    totalCollisionDissipation(step) = totalCollisionDissipation(max(step - 1, 1)) + collisionDissipation;

    % 累加全局碰撞次数
    totalCollisionCount = totalCollisionCount + collisionCount;

    % 实时显示更新
    if mod(step, 100) == 0
        clf;
        scatter3(positions(:, 1), positions(:, 2), positions(:, 3), 20, 'filled');
        xlim([-params.structure_diameter, params.structure_diameter]);
        ylim([-params.structure_diameter, params.structure_diameter]);
        zlim([0, params.structure_height]);
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title(['颗粒模拟 - 当前步数: ', num2str(step)]);
        grid on;
        drawnow;
    end
end

% 打印总碰撞次数
disp(['Total number of collisions across all timesteps: ', num2str(totalCollisionCount)]);
disp('模拟完成！');


% 绘制主结构位移和速度
time = linspace(0, params.numSteps * params.timeStep, params.numSteps);

% 绘制边界碰撞耗能
figure;
plot(time, structureDisplacementHistory, 'b', 'LineWidth', 1.5);
hold on;
xlabel('时间 (s)');
ylabel('主系统振幅 (m)');
title('主系统位移');
legend('主系统位移');
grid on;

% 绘制主系统速度
figure;
plot(time, structureVelocityHistory, 'r', 'LineWidth', 1.5);
hold on;
xlabel('时间 (s)');
ylabel('主系统速度 (m/s)');
title('主系统速度');
legend('主系统速度');
grid on;

% 绘制耗能结果

% 绘制边界碰撞耗能
figure;
plot(time, totalBoundaryDissipation, 'b', 'LineWidth', 1.5);
xlabel('时间 (s)');
ylabel('颗粒-边界碰撞耗能 (J)');
title('颗粒-边界碰撞耗能随时间的变化');
grid on;

% 绘制颗粒间碰撞耗能
figure;
plot(time, totalCollisionDissipation, 'r', 'LineWidth', 1.5);
xlabel('时间 (s)');
ylabel('颗粒-颗粒碰撞耗能 (J)');
title('颗粒-颗粒碰撞耗能随时间的变化');
grid on;

%% 初始化颗粒
function [positions, velocities] = initialize_particles(params)
 % 初始化颗粒的位置和速度
positions = zeros(params.particle_number, 3);

% 圆柱体半径和高度
radius = params.structure_diameter / 2; % 圆柱半径
height = params.structure_height;      % 圆柱高度

for i = 1:params.particle_number
    while true
        % 随机生成 x 和 y 坐标（均匀分布在圆柱体横截面内）
        candidate_x = (rand() - 0.5) * 2 * radius; % x 坐标候选值
        candidate_y = (rand() - 0.5) * 2 * radius; % y 坐标候选值
        
        % 检查是否在圆柱体的径向范围内
        if (candidate_x^2 + candidate_y^2 <= radius^2)
            % 如果在范围内，生成 z 坐标并退出循环
            candidate_z = rand() * height; % z 坐标在 [0, height] 内均匀分布
            positions(i, :) = [candidate_x, candidate_y, candidate_z];
            break;
        end
    end
end


    velocities = zeros(params.particle_number, 3); % 初始速度设置为零
    % velocities = (rand(params.numParticles, 3) - 0.5) * 0.2; % 实验颗粒速度
end

%% 更新颗粒位置
function positions = update_positions(positions, velocities, params)
    % 根据速度更新颗粒位置
    positions = positions + velocities * params.timeStep;
end
%% 应用重力
function velocities = apply_gravity(velocities, params)
    % 为颗粒速度施加重力加速度
    velocities(:, 3) = velocities(:, 3) - params.gravity * params.timeStep;
end
%% 更新主结构动力学
function [structureDisplacement, structureVelocity] = update_structure(externalForce, structureDisplacement, structureVelocity, particlesreactionForce, params)
    % 更新结构的位移和速度（基于牛顿第二定律）
    netForce = externalForce + particlesreactionForce(3); % 合力，包括外力和颗粒的反作用力
    acceleration = (netForce - params.structure_damping * structureVelocity - params.structure_stiffness * structureDisplacement) ...
        / params.structure_mass - params.gravity;

    % 数值积分更新速度和位移
    structureVelocity = structureVelocity + acceleration * params.timeStep;
    structureDisplacement = structureDisplacement + structureVelocity * params.timeStep;
end

%% 颗粒与主结构边界碰撞处理
function [positions, velocities, boundaryDissipation, particlereactionForce] = handle_boundary_collisions(positions, velocities, params, structureVelocity)
    % 初始化能量耗散和反作用力
    boundaryDissipation = 0;
    particlereactionForce = [0, 0, 0];

    % 场景1：底部和顶部碰撞
    bottomCollision = positions(:, 3) < params.particle_radius;
    topCollision = positions(:, 3) > (params.structure_height - params.particle_radius);

    if any(bottomCollision)
        % 底部碰撞处理
        preVelocities = velocities(bottomCollision, 3);
        newVelocities = -params.structure_coefficient_restitution * (preVelocities - structureVelocity) + structureVelocity;
        velocities(bottomCollision, 3) = newVelocities;

        % 冲量计算
        impulse = params.particle_mass * (newVelocities - preVelocities);
        particlereactionForce(3) = particlereactionForce(3) - sum(impulse) / params.timeStep;

        % 能量耗散计算
        preKE = 0.5 * params.particle_mass * preVelocities.^2;
        postKE = 0.5 * params.particle_mass * newVelocities.^2;
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置
        positions(bottomCollision, 3) = params.particle_radius;
    end

    if any(topCollision)
        % 顶部碰撞处理
        preVelocities = velocities(topCollision, 3);
        newVelocities = -params.structure_coefficient_restitution * (preVelocities - structureVelocity) + structureVelocity;
        velocities(topCollision, 3) = newVelocities;

        % 冲量计算
        impulse = params.particle_mass * (newVelocities - preVelocities);
        particlereactionForce(3) = particlereactionForce(3) + sum(impulse) / params.timeStep;

        % 能量耗散计算
        preKE = 0.5 * params.particle_mass * preVelocities.^2;
        postKE = 0.5 * params.particle_mass * newVelocities.^2;
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置
        positions(topCollision, 3) = params.structure_height - params.particle_radius;
    end

    % 场景2：圆柱体壁面碰撞
    cylinderCenter = [0, 0]; % 圆柱中心
    cylinderRadius = params.structure_diameter / 2; % 圆柱半径

    % 计算到圆柱中心的径向距离
    relativePositions = positions(:, 1:2) - cylinderCenter; % 相对圆心位置
    distances = sqrt(sum(relativePositions.^2, 2)); % 径向距离

    % 检测径向碰撞
    radialCollisionMask = distances > (cylinderRadius - params.particle_radius);

    if any(radialCollisionMask)
        % 径向碰撞处理
        preVelocities = velocities(radialCollisionMask, 1:2);

        % 计算径向单位向量
        radialVectors = relativePositions(radialCollisionMask, :) ./ distances(radialCollisionMask);

        % 计算径向速度分量
        radialVelocities = sum(preVelocities .* radialVectors, 2);

        % 碰撞后径向速度
        newRadialVelocities = -params.structure_coefficient_restitution * radialVelocities;

        % 更新速度
        velocityChanges = (newRadialVelocities - radialVelocities) .* radialVectors;
        velocities(radialCollisionMask, 1:2) = velocities(radialCollisionMask, 1:2) + velocityChanges;

        % 冲量计算
        impulse = params.particle_mass * velocityChanges;
        particlereactionForce(1:2) = particlereactionForce(1:2) + sum(impulse, 1) / params.timeStep;

        % 能量耗散计算
        preKE = 0.5 * params.particle_mass * sum(preVelocities.^2, 2);
        postKE = 0.5 * params.particle_mass * sum(velocities(radialCollisionMask, 1:2).^2, 2);
        boundaryDissipation = boundaryDissipation + sum(preKE - postKE);

        % 修正位置
        correctedPositions = cylinderCenter + radialVectors * (cylinderRadius - params.particle_radius);
        positions(radialCollisionMask, 1:2) = correctedPositions;
    end
end

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
        particleVolume = 4/3 * pi * params.particle_diameter^3;
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
                if distance < 2 * params.particle_diameter
                    relVel = velocities(p1, :) - velocities(p2, :);
                    collisionNormal = relPos / max(distance, 1e-6); % 避免零向量
                    relSpeed = dot(relVel, collisionNormal);

                    if relSpeed < 0
                        % 保存碰撞前的速度
                        v1 = velocities(p1, :);
                        v2 = velocities(p2, :);

                        % 计算冲量
                        impulse = -params.particle_coefficient_restitution * relSpeed * collisionNormal;

                        % 更新速度
                        velocities(p1, :) = v1 + 0.5 * impulse;
                        velocities(p2, :) = v2 - 0.5 * impulse;

                        % 碰撞后动能
                        postKE = 0.5 * params.particle_mass * ...
                                 (norm(velocities(p1, :))^2 + norm(velocities(p2, :))^2);

                        % 碰撞后重力势能
                        postPotential = params.particle_mass * g * (positions(p1, 2) + positions(p2, 2));

                        % 总能量
                        preTotalEnergy = 0.5 * params.particle_mass * ...
                                         (norm(v1)^2 + norm(v2)^2) + ...
                                         params.particle_mass * g * (positions(p1, 2) + positions(p2, 2));
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


%% 动态调整时间步长
function timeStep = adjust_time_step(velocities, params)
    % 根据最大速度动态调整时间步长
    maxSpeed = max(sqrt(sum(velocities.^2, 2)));
    timeStep = min(params.timeStep, params.gridSize / maxSpeed * params.safetyFactor);
end
