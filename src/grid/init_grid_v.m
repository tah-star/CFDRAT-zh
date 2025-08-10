function grid_v = init_grid_v(block_info)
%INIT_GRID_V 初始化V速度场的交错网格结构
%   为V动量方程设置交错网格，包括节点分类、边界条件掩码
%   和预组装的常数矩阵系数。V速度位于水平单元面上。
%
%   输入参数:
%       block_info - 障碍物边界信息的元胞数组
%
%   输出参数:
%       grid_v - 包含以下内容的完整网格结构:
%                .is_1st_order/.is_2nd_order - 离散格式掩码
%                .I_const/.J_const/.V_const  - 预组装的边界条件系数
%                .p_idx, .s_idx等            - 节点索引映射
%                (加上几何参数和边界掩码)
%
%   网格布局:
%   - 交错网格: V速度在水平面上(Ny+1 × Nx+2)
%   - 包括左右边界条件的虚拟单元
%   - 上下壁面强制v=0(无穿透条件)
%
%   See also INIT_GRID_U, INIT_GRID_P, GET_AB_V.

global Nx Ny h H L

%% V速度的交错网格几何
% V速度位于压力单元的水平面上
Ny_v = Ny + 1;                       % 每列的水平面数
Nx_v = Nx + 2;                       % 包括左右虚拟列

% V速度节点的物理坐标
xx = h * (-0.5 : 1 : Nx + 0.5);      % 包括虚拟单元: x = -h/2, h/2, ...
yy = h * (0:Ny);                     % 面中心: y = 0, h, 2h, ...
[XX, YY] = meshgrid(xx, yy);

%% 基于几何的节点分类
% 从障碍物几何的主要固体/流体分类
is_solid = get_blocked_mask(block_info, XX, YY);
is_fluid = ~is_solid;

% 浸没边界检测(与固体相邻的流体节点)
is_solid_s = false(Ny_v, Nx_v); is_solid_s(2:end, :) = is_solid(1:end-1, :);
is_solid_n = false(Ny_v, Nx_v); is_solid_n(1:end-1, :) = is_solid(2:end, :);
is_solid_w = false(Ny_v, Nx_v); is_solid_w(:, 2:end) = is_solid(:, 1:end-1);
is_solid_e = false(Ny_v, Nx_v); is_solid_e(:, 1:end-1) = is_solid(:, 2:end);
is_solid_boundary = is_fluid & (is_solid_s | is_solid_n | is_solid_w | is_solid_e);

% 计算域边界识别
is_wall_bottom  = false(Ny_v, Nx_v); is_wall_bottom(1, :) = true;       % 底壁(v=0)
is_wall_top     = false(Ny_v, Nx_v); is_wall_top(end, :) = true;        % 顶壁(v=0)
is_ghost_left   = false(Ny_v, Nx_v); is_ghost_left(2:end-1, 1) = true;  % 左虚拟
is_ghost_right  = false(Ny_v, Nx_v); is_ghost_right(2:end-1, end) = true; % 右虚拟

%% 每个节点的方程类型分配
% 狄利克雷节点: v=0 (壁面和浸没边界的无穿透)
is_dirichlet = is_solid | is_solid_boundary | is_wall_bottom | is_wall_top;

% 诺伊曼虚拟单元: 入口/出口的零梯度
is_neumann = is_ghost_left | is_ghost_right;

% PDE节点: 求解完整动量方程
is_pde_solve = is_fluid & ~is_dirichlet & ~is_neumann;

%% 离散格式选择
% 边界附近使用鲁棒的1阶，内部使用2阶
is_not_pde_solve = ~is_pde_solve;
is_neighbor_s = false(Ny_v, Nx_v); is_neighbor_s(2:end, :) = is_not_pde_solve(1:end-1, :);
is_neighbor_n = false(Ny_v, Nx_v); is_neighbor_n(1:end-1, :) = is_not_pde_solve(2:end, :);
is_neighbor_w = false(Ny_v, Nx_v); is_neighbor_w(:, 2:end) = is_not_pde_solve(:, 1:end-1);
is_neighbor_e = false(Ny_v, Nx_v); is_neighbor_e(:, 1:end-1) = is_not_pde_solve(:, 2:end);
is_near_boundary = is_neighbor_s | is_neighbor_n | is_neighbor_w | is_neighbor_e;

is_1st_order = is_pde_solve & is_near_boundary;   % 边界附近保守
is_2nd_order = is_pde_solve & ~is_near_boundary;  % 内部精确

%% 常数矩阵系数的预组装
% 为稀疏矩阵组装创建全局节点索引
[MM, NN] = meshgrid(1:Nx_v, 1:Ny_v);
p_idx = (NN-1)*Nx_v + MM;                        % 当前节点索引
w_idx = circshift(p_idx, [0, 1]);                % 西邻居
e_idx = circshift(p_idx, [0, -1]);               % 东邻居

% 为稀疏矩阵组装预分配三元组数组
nnz_estimate = nnz(is_dirichlet) + 2*nnz(is_neumann);
I_const = zeros(nnz_estimate, 1);
J_const = zeros(nnz_estimate, 1);
V_const = zeros(nnz_estimate, 1);
current_pos = 0;

% 狄利克雷边界条件: v = 0 (无穿透)
[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_dirichlet, ...
    {p_idx}, 1, ...
    I_const, J_const, V_const, current_pos, p_idx);

% 诺伊曼虚拟单元: ∂v/∂x = 0 → v_ghost = v_interior
[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_left, ...
    {p_idx, e_idx}, [1; -1], ...
    I_const, J_const, V_const, current_pos, p_idx);

[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_right, ...
    {p_idx, w_idx}, [1; -1], ...
    I_const, J_const, V_const, current_pos, p_idx);

% 移除未使用的预分配条目
I_const = I_const(1:current_pos);
J_const = J_const(1:current_pos);
V_const = V_const(1:current_pos);

%% 打包完整的网格结构
grid_v.h = h;
grid_v.L = L;
grid_v.H = H;
grid_v.Nx = Nx_v;
grid_v.Ny = Ny_v;
grid_v.xx = xx;
grid_v.yy = yy;

% 节点分类掩码
grid_v.is_solid = is_solid;
grid_v.is_solid_boundary = is_solid_boundary;
grid_v.is_wall_bottom = is_wall_bottom;
grid_v.is_wall_top = is_wall_top;
grid_v.is_ghost_left = is_ghost_left;
grid_v.is_ghost_right = is_ghost_right;

% 离散格式分类
grid_v.is_1st_order = is_1st_order;
grid_v.is_2nd_order = is_2nd_order;

% 预组装的常数矩阵系数
grid_v.I_const = I_const;
grid_v.J_const = J_const;
grid_v.V_const = V_const;

% 用于模板组装的节点索引映射
grid_v.p_idx = p_idx;
grid_v.s_idx = circshift(p_idx, [1, 0]);
grid_v.n_idx = circshift(p_idx, [-1, 0]);
grid_v.w_idx = w_idx;
grid_v.e_idx = e_idx;
grid_v.ss_idx = circshift(p_idx, [2, 0]);
grid_v.nn_idx = circshift(p_idx, [-2, 0]);
grid_v.ww_idx = circshift(p_idx, [0, 2]);
grid_v.ee_idx = circshift(p_idx, [0, -2]);

end

