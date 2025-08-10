function b = get_b_dp1(u_star, v_star, grid_p)
%GET_B_DP1 组装第一次压力修正方程的右端项
%   通过对中间速度场强制连续性约束，计算PISO算法中第一次压力修正
%   Poisson方程的右端项向量。
%
%   输入参数:
%       u_star, v_star - 预测步骤的中间速度场 [m/s]
%       grid_p         - 压力网格结构
%
%   输出参数:
%       b - 压力修正的右端项向量: ∇²p' = b
%
%   数学背景:
%   动量预测器得到的中间速度u*通常违反连续性方程(∇·u* ≠ 0)。
%   计算压力修正p'将u*投影到散度为零的空间：
%   
%   ∇²p' = (ρ/Δt)∇·u*
%   
%   交错网格上的离散化:
%   - U速度在垂直面上，V速度在水平面上
%   - 散度在单元中心（压力节点）计算
%   - 空间导数使用2阶中心差分
%
%   See also GET_B_DP2, PISO_CORRECT1, PRESSURE_CORRECTION.

global rho dt h

%% 计算压力节点处的速度散度
% 散度计算: ∇·u* = ∂u*/∂x + ∂v*/∂y
% 在交错网格上: div = (u_east - u_west)/h + (v_north - v_south)/h

% X方向速度梯度 (∂u*/∂x)
grad_u_star_x = (u_star(2:end-1, 2:end) - u_star(2:end-1, 1:end-1)) / h;

% Y方向速度梯度 (∂v*/∂y)  
grad_v_star_y = (v_star(2:end, 2:end-1) - v_star(1:end-1, 2:end-1)) / h;

% 每个压力节点的总散度
divergence = grad_u_star_x + grad_v_star_y;

%% 组装右端项向量
% 用(ρ/Δt)缩放散度用于压力修正方程
b_matrix = (rho / dt) * divergence;

% 固体节点的右端项为零（不需要压力修正）
b_matrix(grid_p.is_solid) = 0;

% 负Laplacian矩阵约定的符号修正
b_matrix = -b_matrix;

% 转换为线性求解器的列向量
b_transposed = b_matrix';
b = b_transposed(:);

end
