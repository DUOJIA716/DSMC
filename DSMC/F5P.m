%% 低频段幅频图生成（0-200Hz）
N = length(t);       % 采样点数
Fs = 1/dt;           % 采样频率 (Hz)
freq = Fs*(0:floor(N/2))/N; % 频率轴

% 位移幅频图（0-200Hz）
figure;
hold on;

% 无粒子情况
x_signal = x_no_particles - mean(x_no_particles);
Y = fft(x_signal);
P2 = abs(Y/N);
P1_x_no = P2(1:floor(N/2)+1);
P1_x_no(2:end-1) = 2*P1_x_no(2:end-1);
plot(freq, P1_x_no, 'r-', 'DisplayName', '无颗粒阻尼', 'LineWidth', 1.5);

% DSMC情况
x_signal = x_particles_DSMC - mean(x_particles_DSMC);
Y = fft(x_signal);
P2 = abs(Y/N);
P1_x_DSMC = P2(1:floor(N/2)+1);
P1_x_DSMC(2:end-1) = 2*P1_x_DSMC(2:end-1);
plot(freq, P1_x_DSMC, 'b-', 'DisplayName', '有颗粒阻尼-DSMC', 'LineWidth', 1.5);

% DEM情况
x_signal = x_particles_DEM - mean(x_particles_DEM);
Y = fft(x_signal);
P2 = abs(Y/N);
P1_x_DEM = P2(1:floor(N/2)+1);
P1_x_DEM(2:end-1) = 2*P1_x_DEM(2:end-1);
plot(freq, P1_x_DEM, 'g-', 'DisplayName', '有颗粒阻尼-DEM', 'LineWidth', 1.5);

xlabel('频率 (Hz)');
ylabel('位移幅值 (m)');
title('低频段位移幅频图 (0-200Hz)');
legend;
grid on;
xlim([0 200]); % 聚焦低频段
set(gca, 'FontSize', 10); % 调整字体大小

% 速度幅频图（0-200Hz）
figure;
hold on;

% 无粒子情况
v_signal = v_no_particles - mean(v_no_particles);
Y = fft(v_signal);
P2 = abs(Y/N);
P1_v_no = P2(1:floor(N/2)+1);
P1_v_no(2:end-1) = 2*P1_v_no(2:end-1);
plot(freq, P1_v_no, 'r-', 'DisplayName', '无颗粒阻尼', 'LineWidth', 1.5);

% DSMC情况
v_signal = v_particles_DSMC - mean(v_particles_DSMC);
Y = fft(v_signal);
P2 = abs(Y/N);
P1_v_DSMC = P2(1:floor(N/2)+1);
P1_v_DSMC(2:end-1) = 2*P1_v_DSMC(2:end-1);
plot(freq, P1_v_DSMC, 'b-', 'DisplayName', '有颗粒阻尼-DSMC', 'LineWidth', 1.5);

% DEM情况
v_signal = v_particles_DEM - mean(v_particles_DEM);
Y = fft(v_signal);
P2 = abs(Y/N);
P1_v_DEM = P2(1:floor(N/2)+1);
P1_v_DEM(2:end-1) = 2*P1_v_DEM(2:end-1);
plot(freq, P1_v_DEM, 'g-', 'DisplayName', '有颗粒阻尼-DEM', 'LineWidth', 1.5);

xlabel('频率 (Hz)');
ylabel('速度幅值 (m/s)');
title('低频段速度幅频图 (0-200Hz)');
legend;
grid on;
xlim([0 200]); % 聚焦低频段
set(gca, 'FontSize', 10); % 调整字体大小