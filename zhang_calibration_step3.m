%% 张正友标定法 - 步骤 3: 求解外参 (R, t) 并可视化
clc; clear; close all;

% 1. 加载之前的数据
load('synthetic_calib_data.mat'); % 包含 Ground Truth (extrinsics)
load('homographies.mat');         % 包含 H_all
K_inv = inv(K_true); 

num_images = length(H_all);
estimated_extrinsics = cell(1, num_images);

for i = 1:num_images
    H = H_all{i};
    
    % --- 核心算法: 从 H 分解 R, t ---
    % H = K * [r1, r2, t]  =>  inv(K) * H = [r1, r2, t] (忽略尺度因子)
    h1 = H(:, 1);
    h2 = H(:, 2);
    h3 = H(:, 3);
    
    % 1. 消除尺度因子 lambda
    lambda = 1 / norm(K_inv * h1); 
    
    % 2. 计算 r1, r2, t
    r1 = lambda * K_inv * h1;
    r2 = lambda * K_inv * h2;
    t_est  = lambda * K_inv * h3;
    
    % 3. 计算 r3 (叉乘)
    r3 = cross(r1, r2);
    
    % 4. 组装临时 R
    R_temp = [r1, r2, r3];
    
    % 5. 强制正交化 (SVD方法)
    % 修正噪声，让 R 变成真正的旋转矩阵
    [U, ~, V] = svd(R_temp);
    R_est = U * V';
    
    % 保存结果
    estimated_extrinsics{i}.R = R_est;
    estimated_extrinsics{i}.t = t_est;
    
    % --- 误差分析 ---
    R_true = extrinsics{i}.R;
    t_true = extrinsics{i}.t;
    
    % 计算平移向量的误差 
    t_err = norm(t_est - t_true);
    fprintf('图片 %d: 平移误差 = %.4f mm\n', i, t_err);
end

%% 3D 可视化对比
figure('Color', 'w'); hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('相机位姿恢复结果对比 (红色=真实, 蓝色=估计)');

% 画标定板
plot3(world_points(1,:), world_points(2,:), world_points(3,:), 'k.', 'MarkerSize', 5);
text(0,0,0,'  标定板');

% 画相机
for i = 1:num_images
    % 真实位置 (Ground Truth)
    R_gt = extrinsics{i}.R;
    t_gt = extrinsics{i}.t;
    cam_pos_gt = -R_gt' * t_gt; % 相机在世界坐标系的位置
    
    % 估计位置 (Estimated)
    R_est = estimated_extrinsics{i}.R;
    t_est = estimated_extrinsics{i}.t;
    cam_pos_est = -R_est' * t_est;
    
    % 画图
    plot3(cam_pos_gt(1), cam_pos_gt(2), cam_pos_gt(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    plot3(cam_pos_est(1), cam_pos_est(2), cam_pos_est(3), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
    
    % 连线表示对应关系
    line([cam_pos_gt(1), cam_pos_est(1)], ...
         [cam_pos_gt(2), cam_pos_est(2)], ...
         [cam_pos_gt(3), cam_pos_est(3)], 'Color', 'g');
         
    text(cam_pos_gt(1), cam_pos_gt(2), cam_pos_gt(3), sprintf('  Cam%d', i));
end

view(3); % 设置为3D视角
