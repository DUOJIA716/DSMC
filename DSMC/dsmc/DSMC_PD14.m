% 主程序：基于DSMC方法的颗粒阻尼模拟（含RMS计算、耗能计算和实时显示）
clear; clc;

% 参数设置
params.numParticles = 300; % 颗粒数量
params.boxSize = [0.1, 0.1, 0.4]; % 容器尺寸 [Lx, Ly, Lz]，单位：米
params.particleRadius = 0.00185; % 颗粒半径，单位：米
params.particleMass = 0.01; % 颗粒质量，单位：千克
params.gravity = 9.81; % 重力加速度，单位：m/s^2
params.restitutionCoeff = 0.89; % 恢复系数
params.fillRatioLimit = 0.63; % 最大填充比
params.gridSize = params.particleRadius * 5; % 网格单元边长，单位：米
params.timeStep = 1e-5; % 初始时间步长，单位：秒
params.numSteps = 50000; % 总时间步数
params.safetyFactor = 0.2; % 时间步长调整的安全系数

% 外部荷载参数
params.forceAmplitude = 0.0005 * 1000; % 外部荷载振幅，单位：N
params.forceFrequency = 20; % 外部荷载频率，单位：Hz
params.massStructure = 0.1; % 主体结构质量，单位：kg
params.damping = 0.2; % 主体阻尼系数
params.stiffness = 1000; % 主体刚度，单位：N/m

% 初始化颗粒
[positions, velocities] = initialize_particles(params);
totalEnergyHistory = zeros(params.numSteps, 1); % 保存总能量的历史数据
displacementHistory = zeros(1, params.numSteps); % 记录主系统位移

% 初始化结构运动
structureDisplacement = 0; % 结构初始位移
structureVelocity = 0; % 结构初始速度

% 初始化耗能记录
totalBoundaryDissipation = zeros(1, params.numSteps); % 累积边界碰撞耗能
totalCollisionDissipation = zeros(1, params.numSteps); % 累积颗粒间碰撞耗能


% 初始化全局碰撞次数
totalCollisionCount = 0;

% 模拟循环
for step = 1:params.numSteps
    % 当前时间
    currentTime = step * params.timeStep;

    % 计算外部荷载 F(t) = a * sin(2 * pi * f * t)
    externalForce = params.forceAmplitude * sin(2 * pi * params.forceFrequency * currentTime);

    % 更新结构动力学
    [structureDisplacement, structureVelocity] = update_structure(externalForce, ...
        structureDisplacement, structureVelocity, params);

    % 记录主系统位移
    displacementHistory(step) = structureDisplacement;

    % 自由运动阶段
    positions = update_positions(positions, velocities, params);
    velocities = apply_gravity(velocities, params);

    % 动态调整时间步长
    params.timeStep = adjust_time_step(velocities, params);

    % 边界碰撞处理阶段
    [positions, velocities, boundaryDissipation] = handle_boundary_collisions(positions, velocities, params, structureVelocity);
    totalBoundaryDissipation(step) = totalBoundaryDissipation(max(step - 1, 1)) + boundaryDissipation;

    % 颗粒碰撞处理阶段
    [velocities, collisionDissipation, collisionCount] = handle_collisions(positions, velocities, params);
    totalCollisionDissipation(step) = totalCollisionDissipation(max(step - 1, 1)) + collisionDissipation;

    % 累加全局碰撞次数
    totalCollisionCount = totalCollisionCount + collisionCount;

    % 计算能量并记录
    totalEnergy = calculate_total_energy(positions, velocities, params);
    totalEnergyHistory(step) = totalEnergy;

    % 实时显示更新
    if mod(step, 100) == 0
        clf;
        scatter3(positions(:, 1), positions(:, 2), positions(:, 3), 20, 'filled');
        xlim([0, params.boxSize(1)]);
        ylim([0, params.boxSize(2)]);
        zlim([0, params.boxSize(3)]);
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


% 后处理：计算 RMS 值
rmsValue = sqrt(mean(displacementHistory.^2));
disp(['RMS of primary system amplitude: ', num2str(rmsValue)]);

% 绘制主系统位移和 RMS
time = linspace(0, params.numSteps * params.timeStep, params.numSteps);
figure;
plot(time, displacementHistory, 'b', 'LineWidth', 1.5);
hold on;
yline(rmsValue, '--r', 'LineWidth', 2, 'Label', ['RMS = ', num2str(rmsValue, '%.4f')]);
xlabel('时间 (s)');
ylabel('主系统振幅 (m)');
title('主系统位移和 RMS 值');
legend('主系统位移', 'RMS 值');
grid on;

% 绘制耗能结果
time = linspace(0, params.numSteps * params.timeStep, params.numSteps);

% 绘制边界碰撞耗能
figure;
plot(time, totalBoundaryDissipation, 'b', 'LineWidth', 1.5);
xlabel('时间 (s)');
ylabel('边界碰撞耗能 (J)');
title('边界碰撞耗能随时间的变化');
grid on;

% 绘制颗粒间碰撞耗能
figure;
plot(time, totalCollisionDissipation, 'r', 'LineWidth', 1.5);
xlabel('时间 (s)');
ylabel('颗粒碰撞耗能 (J)');
title('颗粒间碰撞耗能随时间的变化');
grid on;

% 分析和可视化结果
analyze_results(totalEnergyHistory, positions, params);


















