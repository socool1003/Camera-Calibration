clc; clear; close all;

% 1. 相机内参 
width = 640; height = 480;
alpha = 800; % fx
beta = 800;  % fy
u0 = width / 2;
v0 = height / 2;
K_true = [alpha, 0, u0; 0, beta, v0; 0, 0, 1];

fprintf('真实内参矩阵 K:\n');
disp(K_true);

% 2. 棋盘格
% 假设棋盘格在世界坐标系的 Z=0 平面上
square_size = 10; 
board_size = [8, 6]; % 8行6列的角点
[X, Y] = meshgrid(0:board_size(2)-1, 0:board_size(1)-1);
world_points = [X(:)*square_size, Y(:)*square_size, zeros(numel(X), 1)]'; 

% 3. 模拟不同角度拍摄图片 (外参)
num_images = 3;
image_points = cell(1, num_images);
extrinsics = cell(1, num_images);

figure; hold on; grid on;
plot3(world_points(1,:), world_points(2,:), world_points(3,:), 'r.', 'MarkerSize', 10);
title('世界坐标系下的相机位置');
xlabel('X'); ylabel('Y'); zlabel('Z');

for i = 1:num_images
    % 随机生成旋转 (欧拉角) 和 平移
    % 让相机位于标定板上方一定距离，并随机倾斜
    r_angles = (rand(1, 3) - 0.5) * 1.0; % 随机旋转弧度
    R = rotationVectorToMatrix(r_angles); 
    t = [rand*20; rand*20; 1000 + rand*200]; % Z轴距离约1000mm
    
    % 保存外参
    extrinsics{i}.R = R;
    extrinsics{i}.t = t;
    
    % 绘制相机位置
    cam_center = -R' * t;
    plot3(cam_center(1), cam_center(2), cam_center(3), 'bo', 'MarkerSize', 8);
    text(cam_center(1), cam_center(2), cam_center(3), sprintf('Cam %d', i));
    
    % 4. 投影到像素平面 (P = K * [R|t] * P_w)
    P_c = R * world_points + repmat(t, 1, size(world_points, 2)); % 转换到相机坐标系
    P_norm = P_c(1:2, :) ./ P_c(3, :); % 归一化平面坐标 (x/z, y/z)
    
    % 转换到像素坐标
    uv = K_true(1:2, 1:2) * P_norm + K_true(1:2, 3);
    image_points{i} = uv;
end

axis equal;
view(3);

save('synthetic_calib_data.mat', 'world_points', 'image_points', 'K_true', 'extrinsics');

% 罗德里格斯公式
function R = rotationVectorToMatrix(r)
    theta = norm(r);
    if theta < eps
        R = eye(3);
    else
        k = r / theta;
        K = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
        R = eye(3) + sin(theta)*K + (1-cos(theta))*K^2; % 罗德里格斯公式
    end
end