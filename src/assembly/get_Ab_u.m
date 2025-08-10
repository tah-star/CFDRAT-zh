function [A, b] = get_Ab_u(t, u, v, p_pressure, inlet_func, dt, mu, rho, grid_u)
%GET_AB_U 组装U动量方程的系数矩阵和右端项
%   使用混合迎风离散构建U速度分量的线性系统Au = b，
%   并处理时间相关系数。
%
%   输入参数:
%       t           - 当前时间 [s]
%       u, v        - 当前时间步的速度场 [m/s]
%       p_pressure  - 前一时间步的压力场 [Pa]
%       inlet_func  - 入口边界条件函数句柄
%       dt, mu, rho - 时间步长 [s]，动力粘度 [Pa·s]，密度 [kg/m³]
%       grid_u      - 来自init_grid_u的U速度网格结构
%
%   输出参数:
%       A - 动量方程的稀疏系数矩阵
%       b - 包括边界条件和压力梯度项的右端项向量
%
%   离散化方案:
%   - 时间：隐式欧拉法（1阶）
%   - 对流：边界附近1阶迎风，内部2阶
%   - 扩散：2阶中心差分
%   - 压力梯度：2阶中心差分
%
%   See also GET_AB_V, PISO_PREDICT, INTERP_V_TO_U.

%% 提取网格参数和常数
h = grid_u.h;
Ny_u = grid_u.Ny;
Nx_u = grid_u.Nx;

% 节点分类掩码
is_1st_order = grid_u.is_1st_order;
is_2nd_order = grid_u.is_2nd_order;
is_inlet = grid_u.is_inlet;

% 模板组装的节点索引映射
p_idx = grid_u.p_idx; s_idx = grid_u.s_idx; n_idx = grid_u.n_idx;
w_idx = grid_u.w_idx; e_idx = grid_u.e_idx;
ss_idx = grid_u.ss_idx; nn_idx = grid_u.nn_idx;
ww_idx = grid_u.ww_idx; ee_idx = grid_u.ee_idx;

% 预计算常用系数
one_over_dt = 1/dt;
mu_over_h2 = mu/h^2;
one_over_h = 1/h;
one_over_2h = 1/(2*h);

%% 速度插值和迎风因子
% 将V速度插值到U节点位置用于交叉导数项
v_on_u = interp_v_to_u(v);

% 计算迎风离散的速度分量
u_pos = max(u, 0); u_neg = min(u, 0); u_abs = abs(u);
v_on_u_pos = max(v_on_u, 0); v_on_u_neg = min(v_on_u, 0); v_on_u_abs = abs(v_on_u);

%% 有限差分系数矩阵
% 1阶迎风系数（边界附近的鲁棒性）
C1_s = -mu_over_h2 - v_on_u_pos * one_over_h;
C1_w = -mu_over_h2 - u_pos * one_over_h;
C1_p = one_over_dt + u_abs * one_over_h + v_on_u_abs * one_over_h + 4*mu_over_h2;
C1_e = -mu_over_h2 + u_neg * one_over_h;
C1_n = -mu_over_h2 + v_on_u_neg * one_over_h;

% 2阶迎风系数（内部精确）
C2_ss = v_on_u_pos * one_over_2h;
C2_s  = -mu_over_h2 - 2*v_on_u_pos * one_over_h;
C2_ww = u_pos * one_over_2h;
C2_w  = -mu_over_h2 - 2*u_pos * one_over_h;
C2_p  = one_over_dt + 3*u_abs * one_over_2h + 3*v_on_u_abs * one_over_2h + 4*mu_over_h2;
C2_e  = -mu_over_h2 + 2*u_neg * one_over_h;
C2_ee = -u_neg * one_over_2h;
C2_n  = -mu_over_h2 + 2*v_on_u_neg * one_over_h;
C2_nn = -v_on_u_neg * one_over_2h;

%% 使用IJV三元组格式的稀疏矩阵组装
% 估算非零元素数量以进行高效内存分配
nnz_const = length(grid_u.I_const);
nnz_1st = 5 * nnz(is_1st_order);
nnz_2nd = 9 * nnz(is_2nd_order);
total_nnz = nnz_const + nnz_1st + nnz_2nd;

I = zeros(total_nnz, 1);
J = zeros(total_nnz, 1);
V = zeros(total_nnz, 1);

% 插入预计算的常数边界条件系数
I(1:nnz_const) = grid_u.I_const;
J(1:nnz_const) = grid_u.J_const;
V(1:nnz_const) = grid_u.V_const;
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
NyNx = Ny_u * Nx_u;
A = sparse(I, J, V, NyNx, NyNx);

%% 右端项向量组装
b_matrix = zeros(Ny_u, Nx_u);

% 入口边界条件：时间相关的抛物线剖面
yy_inlet = grid_u.yy(is_inlet);
u_inlet_values = inlet_func(t, yy_inlet);
b_matrix(is_inlet) = u_inlet_values;
% 其他边界条件节点默认为零

% 内部节点：时间项减去压力梯度
grad_p_x = (p_pressure(:, 2:end) - p_pressure(:, 1:end-1)) / h;
u_internal = u(2:end-1, 2:end-1);
b_matrix(2:end-1, 2:end-1) = u_internal * one_over_dt - grad_p_x / rho;

% 强制固体节点的右端项为零
b_matrix(grid_u.is_solid) = 0;
b_matrix(grid_u.is_solid_boundary) = 0;

% 转换为列向量（MATLAB稀疏矩阵约定）
b_transposed = b_matrix';
b = b_transposed(:);

end


function v_on_u = interp_v_to_u(v)
%INTERP_V_TO_U 将V速度从V节点插值到U节点位置
%   对内部节点使用4点平均，对虚拟单元进行外推，
%   假设无滑移边界条件。

    % 内部U节点的4点平均
    v_on_u_interior = (v(1:end-1, 1:end-1) + v(1:end-1, 2:end) + ...
                       v(2:end, 1:end-1)   + v(2:end, 2:end)) / 4;
    
    % 来自无滑移条件的虚拟单元值：(v_ghost + v_interior)/2 = 0
    v_on_u_ghost_bottom = -v_on_u_interior(1, :);
    v_on_u_ghost_top    = -v_on_u_interior(end, :);
    
    v_on_u = [v_on_u_ghost_bottom; v_on_u_interior; v_on_u_ghost_top];
end
