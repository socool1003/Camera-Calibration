%% 张正友标定法 - 步骤 4: 非线性优化 (Bundle Adjustment)
clc; clear; close all;

% 1. 加载所有初值 
load('synthetic_calib_data.mat');      % 真实值
load('homographies.mat');              % H 矩阵

K_init = [800.0000, 0, 320.0000; 0, 800.0000, 240.0000; 0, 0, 1]; 

K_init(1,1) = K_init(1,1) + 5; 
K_init(2,2) = K_init(2,2) + 5;

num_images = length(extrinsics);
r_vecs_init = zeros(3, num_images);
t_vecs_init = zeros(3, num_images);
for i = 1:num_images
    % 使用真实外参作为初值 (使用 Step 3 的结果)
    R = extrinsics{i}.R;
    t = extrinsics{i}.t;
    r_vecs_init(:, i) = rotationMatrixToVector(R); % 转为旋转向量
    t_vecs_init(:, i) = t;
end
dist_coeffs_init = [0; 0]; % 初始假设无畸变 [k1, k2]
% -----------------------------------------------------------

fprintf('开始非线性优化 (Bundle Adjustment)...\n');

% 2. 参数打包 (Vectorization)
% lsqnonlin 需要把所有待优化参数把它拉成一个长向量
% 参数结构: [fx, fy, u0, v0, k1, k2, r1_x, r1_y, r1_z, t1_x, ..., rn_..., tn_...]
x0 = pack_params(K_init, dist_coeffs_init, r_vecs_init, t_vecs_init);

% 3. 定义优化选项
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', ...
                       'Display', 'iter', 'MaxFunctionEvaluations', 10000);

% 4. 运行优化
% world_points: 3xN, image_points: cell array
[x_opt, resnorm, residual, exitflag] = lsqnonlin(@(x) cost_function(x, world_points, image_points), x0, [], [], options);

% 5. 解包参数
[K_opt, dist_opt, r_vecs_opt, t_vecs_opt] = unpack_params(x_opt, num_images);

% 结果展示与对比
fprintf('\n--- 优化结果对比 ---\n');
fprintf('真实 K: fx=%.2f, fy=%.2f\n', K_true(1,1), K_true(2,2));
fprintf('初始 K: fx=%.2f, fy=%.2f\n', K_init(1,1), K_init(2,2));
fprintf('优化 K: fx=%.2f, fy=%.2f\n', K_opt(1,1), K_opt(2,2));

fprintf('\n畸变系数 [k1, k2]:\n');
fprintf('优化后: [%.5f, %.5f]\n', dist_opt(1), dist_opt(2));

% 计算重投影误差
initial_error = compute_reprojection_error(x0, world_points, image_points);
final_error = compute_reprojection_error(x_opt, world_points, image_points);

fprintf('\n平均重投影误差 (像素):\n');
fprintf('优化前: %.4f pixels\n', initial_error);
fprintf('优化后: %.4f pixels\n', final_error);

% 绘图: 重投影误差可视化
figure('Color', 'w');
subplot(1,2,1);
bar([initial_error, final_error]);
set(gca, 'XTickLabel', {'优化前', '优化后'});
title('平均重投影误差对比');
ylabel('像素误差');

subplot(1,2,2);
% 画出最后一张图的重投影点
img_idx = 1;
projected_pts = project_points(world_points, K_opt, dist_opt, r_vecs_opt(:,img_idx), t_vecs_opt(:,img_idx));
measured_pts = image_points{img_idx};

plot(measured_pts(1,:), measured_pts(2,:), 'go', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
plot(projected_pts(1,:), projected_pts(2,:), 'r+', 'MarkerSize', 8, 'LineWidth', 1.5);
legend('真实观测点', '优化后投影点');
title(['第1张图的重投影效果']);
axis equal; grid on; set(gca, 'YDir', 'reverse');

function residuals = cost_function(x, world_points, image_points)
    num_images = length(image_points);
    [K, dist, r_vecs, t_vecs] = unpack_params(x, num_images);
    
    residuals = [];
    for i = 1:num_images
        % 计算投影
        proj_uv = project_points(world_points, K, dist, r_vecs(:,i), t_vecs(:,i));
        % 观测值
        meas_uv = image_points{i};
        % 误差 (2xN -> 向量)
        diff = proj_uv - meas_uv;
        residuals = [residuals; diff(:)];
    end
end

function uv = project_points(P_world, K, dist, r_vec, t_vec)
    % 1. 世界 -> 相机
    R = rotationVectorToMatrix_local(r_vec);
    P_cam = R * P_world + t_vec;
    
    % 2. 归一化平面
    x = P_cam(1, :) ./ P_cam(3, :);
    y = P_cam(2, :) ./ P_cam(3, :);
    
    % 3. 畸变 (Radial Distortion)
    k1 = dist(1); k2 = dist(2);
    r2 = x.^2 + y.^2;
    dist_factor = (1 + k1*r2 + k2*r2.^2);
    x_dist = x .* dist_factor;
    y_dist = y .* dist_factor;
    
    % 4. 像素坐标
    u = K(1,1) * x_dist + K(1,3);
    v = K(2,2) * y_dist + K(2,3);
    uv = [u; v];
end

% 参数打包
function x = pack_params(K, dist, r_vecs, t_vecs)
    x = [K(1,1); K(2,2); K(1,3); K(2,3); dist; r_vecs(:); t_vecs(:)];
end

% 参数解包
function [K, dist, r_vecs, t_vecs] = unpack_params(x, num_images)
    fx = x(1); fy = x(2); u0 = x(3); v0 = x(4);
    K = [fx, 0, u0; 0, fy, v0; 0, 0, 1];
    dist = x(5:6);
    
    offset = 6;
    num_ext = num_images * 3;
    r_vecs = reshape(x(offset+1 : offset+num_ext), 3, num_images);
    t_vecs = reshape(x(offset+num_ext+1 : end), 3, num_images);
end

function mean_err = compute_reprojection_error(x, world_points, image_points)
    res = cost_function(x, world_points, image_points);
    mean_err = mean(abs(res)); % 简单的平均绝对误差
end

% 罗德里格斯公式
function R = rotationVectorToMatrix_local(r)
    theta = norm(r);
    if theta < eps
        R = eye(3);
    else
        k = r / theta;
        K_cross = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
        R = eye(3) + sin(theta)*K_cross + (1-cos(theta))*K_cross^2;
    end
end

function r = rotationMatrixToVector(R)
    % 简单逆变换 (trace相关)
    theta = acos((trace(R)-1)/2);
    if theta < eps
        r = [0;0;0];
    else
        r = theta * 1/(2*sin(theta)) * [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
    end
end