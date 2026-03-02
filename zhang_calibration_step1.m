%% 张正友标定法 - 步骤 1: 计算单应性矩阵 H
clc; clear; 
load('synthetic_calib_data.mat'); 

num_images = length(image_points);
H_all = cell(1, num_images); % 存放每一张图片的 H 矩阵

% 取世界坐标的 X, Y (Z=0)
M = world_points(1:2, :); % 2xN 矩阵
num_points = size(M, 2);

for i = 1:num_images
    m = image_points{i}; % 当前图片的像素坐标 2xN
    
    % 构建 DLT (直接线性变换) 方程组 A*h = 0
    A = zeros(2 * num_points, 9);
    
    for j = 1:num_points
        X = M(1, j); Y = M(2, j);
        u = m(1, j); v = m(2, j);
        
        % 每一对点提供两个方程
        A(2*j-1, :) = [X, Y, 1, 0, 0, 0, -u*X, -u*Y, -u];
        A(2*j, :)   = [0, 0, 0, X, Y, 1, -v*X, -v*Y, -v];
    end
    
    % SVD 求解 A*h = 0 的最小二乘解
    % 解是 A 的最小奇异值对应的右奇异向量 (V 的最后一列)
    [~, ~, V] = svd(A);
    h = V(:, end);
    
    % 重组为 3x3 矩阵
    H = reshape(h, 3, 3)';
    
    % 归一化 
    H = H / H(3, 3);
    
    H_all{i} = H;
end

save('homographies.mat', 'H_all');
disp('所有 H 矩阵已保存。');