function [A, b] = get_Ab_v(u, v, p_pressure, dt, mu, rho, grid_v)
%GET_AB_V 组装V动量方程的系数矩阵和右端项
%   使用混合迎风离散构建V速度分量的线性系统Av = b，
%   并处理时间相关系数。
%
%   输入参数:
%       u, v        - 当前时间步的速度场 [m/s]
%       p_pressure  - 前一时间步的压力场 [Pa]
%       dt, mu, rho - 时间步长 [s]，动力粘度 [Pa·s]，密度 [kg/m³]
%       grid_v      - 来自init_grid_v的V速度网格结构
%
%   输出参数:
%       A - 动量方程的稀疏系数矩阵
%       b - 包括压力梯度项的右端项向量
%
%   离散化方案:
%   - 时间：隐式欧拉法（1阶）
%   - 对流：边界附近1阶迎风，内部2阶
%   - 扩散：2阶中心差分
%   - 压力梯度：2阶中心差分（y方向）
%
%   注意：V速度在上下壁面有无穿透边界条件（v=0），
%   在入口/出口有零梯度边界条件。
%
%   See also GET_AB_U, PISO_PREDICT, INTERP_U_TO_V.

%% 提取网格参数和常数
h = grid_v.h;
Ny_v = grid_v.Ny;
Nx_v = grid_v.Nx;

% 节点分类掩码
is_1st_order = grid_v.is_1st_order;
is_2nd_order = grid_v.is_2nd_order;

% 模板组装的节点索引映射
p_idx = grid_v.p_idx; s_idx = grid_v.s_idx; n_idx = grid_v.n_idx;
w_idx = grid_v.w_idx; e_idx = grid_v.e_idx;
ss_idx = grid_v.ss_idx; nn_idx = grid_v.nn_idx;
ww_idx = grid_v.ww_idx; ee_idx = grid_v.ee_idx;

% 预计算常用系数
one_over_dt = 1/dt;
mu_over_h2 = mu/h^2;
one_over_h = 1/h;
one_over_2h = 1/(2*h);

%% 速度插值和迎风因子
% 将U速度插值到V节点位置用于交叉导数项
u_on_v = interp_u_to_v(u);

% 计算迎风离散的速度分量
u_on_v_pos = max(u_on_v, 0); u_on_v_neg = min(u_on_v, 0); u_on_v_abs = abs(u_on_v);
v_pos = max(v, 0); v_neg = min(v, 0); v_abs = abs(v);

%% 有限差分系数矩阵
% 1阶迎风系数（边界附近的鲁棒性）
C1_s = -mu_over_h2 - v_pos * one_over_h;
C1_w = -mu_over_h2 - u_on_v_pos * one_over_h;
C1_p = one_over_dt + u_on_v_abs * one_over_h + v_abs * one_over_h + 4*mu_over_h2;
C1_e = -mu_over_h2 + u_on_v_neg * one_over_h;
C1_n = -mu_over_h2 + v_neg * one_over_h;

% 2阶迎风系数（内部精确）
C2_ss = v_pos * one_over_2h;
C2_s  = -mu_over_h2 - 2*v_pos * one_over_h;
C2_ww = u_on_v_pos * one_over_2h;
C2_w  = -mu_over_h2 - 2*u_on_v_pos * one_over_h;
C2_p  = one_over_dt + 3*u_on_v_abs * one_over_2h + 3*v_abs * one_over_2h + 4*mu_over_h2;
C2_e  = -mu_over_h2 + 2*u_on_v_neg * one_over_h;
C2_ee = -u_on_v_neg * one_over_2h;
C2_n  = -mu_over_h2 + 2*v_neg * one_over_h;
C2_nn = -v_neg * one_over_2h;

%% 使用IJV三元组格式的稀疏矩阵组装
% 估算非零元素数量以进行高效内存分配
nnz_const = length(grid_v.I_const);
nnz_1st = 5 * nnz(is_1st_order);
nnz_2nd = 9 * nnz(is_2nd_order);
total_nnz = nnz_const + nnz_1st + nnz_2nd;

I = zeros(total_nnz, 1);
J = zeros(total_nnz, 1);
V = zeros(total_nnz, 1);

% 插入预计算的常数边界条件系数
I(1:nnz_const) = grid_v.I_const;
J(1:nnz_const) = grid_v.J_const;
V(1:nnz_const) = grid_v.V_const;
current_pos = nnz_const;

% 为边界邻近节点组装1阶模板
[I, J, V, current_pos] = assemble_stencil_dynamic(is_1st_order, ...
    {s_idx, w_idx, p_idx, e_idx, n_idx}, ...
    {C1_s, C1_w, C1_p, C1_e, C1_n}, ...
    I, J, V, current_pos, p_idx);

% 为内部节点组装2阶模板
[I, J, V, current_pos] = assemble_stencil_dynamic(is_2nd_order, ...
    {ss_idx, s_idx, ww_idx, w_idx, p_idx, e_idx, ee_idx, n_idx, nn_idx}, ...
    {C2_ss, C2_s, C2_ww, C2_w, C2_p, C2_e, C2_ee, C2_n, C2_nn}, ...
    I, J, V, current_pos, p_idx);

% 创建最终稀疏矩阵
NyNx = Ny_v * Nx_v;
A = sparse(I, J, V, NyNx, NyNx);

%% 右端项向量组装
b_matrix = zeros(Ny_v, Nx_v);
% V速度的所有边界条件默认为零

% 内部节点：时间项减去压力梯度（y方向）
grad_p_y = (p_pressure(2:end, :) - p_pressure(1:end-1, :)) / h;
v_internal = v(2:end-1, 2:end-1);
b_matrix(2:end-1, 2:end-1) = v_internal * one_over_dt - grad_p_y / rho;

% 强制固体节点的右端项为零
b_matrix(grid_v.is_solid) = 0;
b_matrix(grid_v.is_solid_boundary) = 0;

% 转换为列向量（MATLAB稀疏矩阵约定）
b_transposed = b_matrix';
b = b_transposed(:);

end

function u_on_v = interp_u_to_v(u)
%INTERP_U_TO_V 将U速度从U节点插值到V节点位置
%   对内部节点使用4点平均，适当处理入口（Dirichlet）
%   和出口（Neumann）条件的边界。

    % 内部V节点的4点平均
    u_on_v_interior = (u(1:end-1, 1:end-1) + u(1:end-1, 2:end) + ...
                       u(2:end, 1:end-1)   + u(2:end, 2:end)) / 4;
    
    % 出口：零梯度条件（Neumann）
    u_on_v_ghost_right = u_on_v_interior(:, end);

    % 入口：从Dirichlet条件外推
    u_left_bound = (u(2:end,1) + u(1:end-1,1)) / 2; 
    u_on_v_ghost_left  =  2 * u_left_bound - u_on_v_interior(:, 1);

    u_on_v = [u_on_v_ghost_left, u_on_v_interior, u_on_v_ghost_right];
end

