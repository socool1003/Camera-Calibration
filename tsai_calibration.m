%% 蔡氏标定法 (Tsai's Method)
clc; clear; close all;

% 1. 加载数据
load('synthetic_calib_data.mat');

% 选取第 1 张图片的数据进行测试
img_idx = 1;
P_world = world_points;       % 3xN
p_img   = image_points{img_idx}; % 2xN

% 真实值 
R_true = extrinsics{img_idx}.R;
t_true = extrinsics{img_idx}.t;
f_true = K_true(1,1); % fx

% 假设主点已知
u0 = K_true(1,3);
v0 = K_true(2,3);

% 图像坐标转换为去中心的图像物理坐标
% x_d = u - u0; y_d = v - v0;
xd = p_img(1, :) - u0;
yd = p_img(2, :) - v0;

num_points = length(xd);

%第一步: 利用 RAC 求解 R 和 tx, ty
% RAC (Radial Alignment Constraint): 
% 径向排列约束原理：图像点到主点的向量方向，不受径向畸变影响。
% 公式: x_d / y_d = x_u / y_u = (r1*X + r2*Y + tx) / (r4*X + r5*Y + ty)
% 整理得线性方程: [v*X, v*Y, -u*X, -u*Y, v] * [r1/ty; r2/ty; r4/ty; r5/ty; tx/ty] = u

% 构建线性方程组 A * L = b
A = [];
b = [];

for i = 1:num_points
    X = P_world(1, i);
    Y = P_world(2, i);
    u = xd(i);
    v = yd(i);
    
    % 为了数值稳定性，剔除 v 接近 0 的点
    if abs(v) > 1e-5
        row = [v*X, v*Y, -u*X, -u*Y, v];
        A = [A; row];
        b = [b; u];
    end
end

% 最小二乘求解 L
L = A \ b;

% L = [r1/ty; r2/ty; r4/ty; r5/ty; tx/ty]
% 从 L 中恢复参数
ty_abs = 1 / sqrt(L(1)^2 + L(2)^2 + L(3)^2 + L(4)^2); % 此时 ty 只有大小，符号未知

% 尝试 ty 的正负号 (通常根据物体在相机前方确定)
% 这里简化处理，选取使结果合理的符号
scale = ty_abs; 

% 恢复出的中间变量
r1 = L(1) * scale;
r2 = L(2) * scale;
r4 = L(3) * scale;
r5 = L(4) * scale;
tx = L(5) * scale;
ty = scale; % 暂时假设 ty > 0

% 组装旋转矩阵的前两列
R_col1 = [r1; r4; sqrt(1 - r1^2 - r4^2)]; % 利用单位向量约束求 r7
R_col2 = [r2; r5; sqrt(1 - r2^2 - r5^2)]; % 利用单位向量约束求 r8

% 这里的 R 求解非常粗糙，Tsai 原文有更复杂的 SVD 修正
% 为了演示，我们直接用正交化修正一下
R_temp = [R_col1, R_col2, cross(R_col1, R_col2)];
[U_r, ~, V_r] = svd(R_temp);
R_est = U_r * V_r'; % 强制正交化

t_est_xy = [tx; ty]; 

fprintf('--- 第一步结果 (RAC) ---\n');
fprintf('真实 R(1:2,1:2):\n'); disp(R_true(1:2,1:2));
fprintf('估计 R(1:2,1:2):\n'); disp(R_est(1:2,1:2));

%第二步: 求解 tz 和 f (焦距) 
% 此时 R, tx, ty 已知。
% 公式: y_u = f * Y_c / Z_c
% y_d = f * (r4*X + r5*Y + ty) / (r7*X + r8*Y + tz)
% 整理为线性方程: [y_d, -(r4*X + r5*Y + ty)] * [tz; f] = y_d * (r7*X + r8*Y)

A2 = [];
b2 = [];

r4 = R_est(2,1); r5 = R_est(2,2);
r7 = R_est(3,1); r8 = R_est(3,2);

for i = 1:num_points
    X = P_world(1, i);
    Y = P_world(2, i);
    v = yd(i);
    
    Yc_partial = r4*X + r5*Y + ty;
    Zc_partial = r7*X + r8*Y;
    
    A2 = [A2; v, -Yc_partial];
    b2 = [b2; v * Zc_partial];
end

% 求解 [tz; f]
res = A2 \ b2;
tz_est = res(1);
f_est = res(2);

t_est = [tx; ty; tz_est];

% 结果对比
fprintf('\n--- 最终结果对比 (Tsai) ---\n');
fprintf('焦距 f:\n  真实值: %.2f\n  估计值: %.2f\n', f_true, f_est);
fprintf('平移向量 t:\n');
disp([t_true, t_est]);
fprintf('平移误差: %.4f mm\n', norm(t_true - t_est));

% 3D 可视化 
figure('Color', 'w'); hold on; grid on; axis equal;
plot3(P_world(1,:), P_world(2,:), P_world(3,:), 'k.');
text(0,0,0, '标定板');

% 画真实相机
C_true = -R_true' * t_true;
plot3(C_true(1), C_true(2), C_true(3), 'go', 'MarkerSize', 10, 'LineWidth', 2);

% 画 Tsai 估计相机
C_est = -R_est' * t_est;
plot3(C_est(1), C_est(2), C_est(3), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

legend('世界坐标点', '真实位置', 'Tsai估计位置');
title('Tsai 两步法标定结果');
view(3);
