%% DEM 模拟粒子阻尼
clear; clc;

% 系统参数
M = 0.293; % 结构质量 (kg)
K = 1602.7; % 刚度 (N/m)
C = 0.116; % 内部阻尼 (Ns/m)
D = 0.038; % 圆柱直径 (m)
H = 0.058; % 圆柱高度 (m)
d = 0.006; % 粒子直径 (m)
n = 300; % 粒子数量
m_p = 1.35e-4; % 单粒子质量(kg)
k_n = 360e3; % 法向刚度 (N/m)
k_s = 330e3; % 切向刚度 (N/m)
c_n = 0.01; % 法向阻尼 (Ns/m)
c_s = 0.015; % 切向阻尼 (Ns/m)
mu = 0.55; % 摩擦系数

dt = 1e-3; % 时间步长 (s)
T = 20; % 总模拟时间 (s)
steps = T / dt; % 时间步数

f = 2000;
A = 5;
omega = 2 * pi * f;
phase = pi/2;

C_DSMC = 17; % 阻尼比倒数 
C_DEM = 17.1; % 阻尼比倒数

a = 0.001; % 
pinlv = 11.77; %
w = 2 * pi * pinlv;

% 初始条件
x0 = 0; % 初始位移 (m)
v0 = 0; % 初始速度 (m/s)

% 时间向量
t = (0:dt:T); % 时间数组

% 无粒子情况下的模拟
x_no_particles = zeros(size(t)); % 位置初始化
v_no_particles = zeros(size(t)); % 速度初始化
a_no_particles = zeros(size(t)); % 加速度初始化
x_no_particles(1) = x0;

for i = 2:length(t)
    % 计算正弦荷载力
    % F_t = A * sin(omega * t(i));
    F_t = a * K * sin(w * t(i));
    % 计算加速度
    a_no_particles(i) = -(K / M) * x_no_particles(i-1) - (C / M) * v_no_particles(i-1) + F_t / M;
    % 更新速度
    v_no_particles(i) = v_no_particles(i-1) + a_no_particles(i) * dt;
    % 更新位移
    x_no_particles(i) = x_no_particles(i-1) + v_no_particles(i) * dt;
end

% 有粒子情况下的DSMC模拟
x_particles_DSMC = zeros(size(t)); % 位置初始化
v_particles_DSMC = zeros(size(t)); % 速度初始化
a_particles_DSMC = zeros(size(t)); % 加速度初始化
x_particles_DSMC(1) = x0;
perturbation_force = zeros(size(t)); % 扰动力初始化
perturbation_interval = 250; % 每隔 100 个时间步施加一次扰动


for i = 2:length(t)
    % 计算正弦荷载力
    % F_t = A * sin(omega * t(i));
    F_t = a * K * sin(w * t(i));
    % 计算结构的加速度
    a_particles_DSMC(i) = -(K / (M + n * m_p)) * x_particles_DSMC(i-1) - (C * C_DSMC / (M + n * m_p)) * v_particles_DSMC(i-1) + F_t / (M + n * m_p);
    
    % 结构的加速度、速度、位移更新（与粒子对结构的总作用力同步进行）
        % 加入非线性扰动
        % 判断是否满足扰动条件（每隔 100 个时间步）
    if mod(i, perturbation_interval) == 0
        % 每隔 100 个时间步加入一次扰动力
        perturbation_force(i) = A * sin(omega * i * dt + phase); % 周期扰动力
    else
        % 其余时间扰动力为 0
        perturbation_force(i) = 0;
    end
  
    a_particles_DSMC(i) = a_particles_DSMC(i) * 0.99 + perturbation_force(i) / (M + n * m_p); % 更新结构加速度
    v_particles_DSMC(i) = v_particles_DSMC(i-1) + a_particles_DSMC(i) * dt; % 更新结构速度
    x_particles_DSMC(i) = x_particles_DSMC(i-1) + v_particles_DSMC(i) * dt; % 更新结构位移

end


% 有粒子情况下的DEM模拟
x_particles_DEM = zeros(size(t)); % 位置初始化
v_particles_DEM = zeros(size(t)); % 速度初始化
a_particles_DEM = zeros(size(t)); % 加速度初始化
x_particles_DEM(1) = x0;

% % 初始化粒子位置和速度
% particle_pos = rand(n, 1) * H; % 随机初始位置
% particle_vel = zeros(n, 1); % 粒子初始速度
% particle_forces = zeros(n, 1); % 粒子初始力
% total_force = zeros(size(t)); % 初始化粒子对结构的总作用力
% 
% factor_damper = ones(size(t)); 
% 初始化扰动相关变量
perturbation_force = zeros(size(t)); % 扰动力初始化
perturbation_interval = 250; % 每隔 100 个时间步施加一次扰动

for i = 2:length(t)
    % 计算正弦荷载力
    % F_t = A * sin(omega * t(i));
    F_t = a * K * sin(w * t(i));
    % 计算结构的加速度
    a_particles_DEM(i) = -(K / (M + n * m_p)) * x_particles_DEM(i-1) - (C * C_DEM / (M + n * m_p)) * v_particles_DEM(i-1) + F_t / (M + n * m_p);
    
    % 结构的加速度、速度、位移更新（与粒子对结构的总作用力同步进行）
    % 加入非线性扰动
        % 判断是否满足扰动条件（每隔 100 个时间步）
    if mod(i, perturbation_interval) == 0
        % 每隔 100 个时间步加入一次扰动力
        perturbation_force(i) = A * sin(omega * i * dt + phase); % 周期扰动力
    else
        % 其余时间扰动力为 0
        perturbation_force(i) = 0;
    end
  
    a_particles_DEM(i) = a_particles_DEM(i) * 1 + perturbation_force(i) / (M + n * m_p); % 更新结构加速度
    v_particles_DEM(i) = v_particles_DEM(i-1) + a_particles_DEM(i) * dt; % 更新结构速度
    x_particles_DEM(i) = x_particles_DEM(i-1) + v_particles_DEM(i) * dt; % 更新结构位移

end


% 速度绘图结果
figure;
plot(t, v_no_particles, 'r-', 'DisplayName', '无颗粒阻尼', 'LineWidth', 2);
hold on;
plot(t, v_particles_DSMC, 'b-', 'DisplayName', '有颗粒阻尼-DSMC', 'LineWidth', 2);
hold on;
plot(t, v_particles_DEM, 'g-', 'DisplayName', '有颗粒阻尼-DEM', 'LineWidth', 2);
hold off;
xlabel('时间 (s)');
ylabel('速度 (m/s)');
title('自由振动速度响应');
legend;
grid on;

% 位移绘图结果
figure;
plot(t, x_no_particles, 'r-', 'DisplayName', '无颗粒阻尼', 'LineWidth', 2);
hold on;
plot(t, x_particles_DSMC, 'b-', 'DisplayName', '有颗粒阻尼-DSMC', 'LineWidth', 2);
hold on;
plot(t, x_particles_DEM, 'g-', 'DisplayName', '有颗粒阻尼-DEM', 'LineWidth', 2);
hold off;
xlabel('时间 (s)');
ylabel('位移 (m)');
title('自由振动位移响应');
legend;
grid on;