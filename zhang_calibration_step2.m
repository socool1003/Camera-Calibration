%% 张正友标定法 - 步骤 2: 求解内参 K
clc; clear;
load('homographies.mat');

% 定义 Vij 函数
% v_ij 表示 H 的第 i 列和第 j 列的特定组合
get_v = @(H, i, j) [H(1,i)*H(1,j), ...
                    H(1,i)*H(2,j) + H(2,i)*H(1,j), ...
                    H(2,i)*H(2,j), ...
                    H(3,i)*H(1,j) + H(1,i)*H(3,j), ...
                    H(3,i)*H(2,j) + H(2,i)*H(3,j), ...
                    H(3,i)*H(3,j)]';

% 构建方程组 V*b = 0
% b 是一个 6 维向量，包含了 K 的参数信息
% 每一张图片提供 2 个约束方程
num_images = length(H_all);
V_mat = [];

for i = 1:num_images
    H = H_all{i};
    % 约束 1: h1 和 h2 正交 -> h1^T * B * h2 = 0
    v12 = get_v(H, 1, 2);
    
    % 约束 2: h1 和 h2 模长相等 -> h1^T * B * h1 - h2^T * B * h2 = 0
    v11 = get_v(H, 1, 1);
    v22 = get_v(H, 2, 2);
    
    V_mat = [V_mat; v12'; (v11 - v22)'];
end

% SVD 求解 V_mat * b = 0
[~, ~, V_svd] = svd(V_mat);
b = V_svd(:, end);

% 从 b 中提取参数 (Closed-form solution)
% b = [B11, B12, B22, B13, B23, B33]
B11 = b(1); B12 = b(2); B22 = b(3); B13 = b(4); B23 = b(5); B33 = b(6);

% 根据张正友论文附录中的公式，反解内参
v0 = (B12*B13 - B11*B23) / (B11*B22 - B12^2);
lambda = B33 - (B13^2 + v0*(B12*B13 - B11*B23)) / B11;
alpha = sqrt(lambda / B11);
beta = sqrt(lambda * B11 / (B11*B22 - B12^2));
gamma = -B12 * alpha^2 * beta / lambda;
u0 = gamma * v0 / beta - B13 * alpha^2 / lambda;

% 组装计算出的 K
K_estimated = [alpha, gamma, u0; 
               0,     beta,  v0; 
               0,     0,     1];

fprintf('\n--- 标定结果 ---\n');
disp('计算出的内参矩阵 K:');
disp(K_estimated);

load('synthetic_calib_data.mat', 'K_true');
disp('真实内参矩阵 K (Ground Truth):');
disp(K_true);

% 计算误差
error_matrix = abs(K_estimated - K_true);
fprintf('绝对误差:\n');
disp(error_matrix);