function grid_u = init_grid_u(block_info)
%INIT_GRID_U 初始化U速度场的交错网格结构
%   为U动量方程设置交错网格，包括节点分类、边界条件掩码
%   和预组装的常数矩阵系数。U速度位于垂直单元面上。
%
%   输入参数:
%       block_info - 障碍物边界信息的元胞数组
%
%   输出参数:
%       grid_u - 包含以下内容的完整网格结构:
%                .is_1st_order/.is_2nd_order - 离散格式掩码
%                .I_const/.J_const/.V_const  - 预组装的边界条件系数
%                .p_idx, .s_idx等            - 节点索引映射
%                (加上几何参数和边界掩码)
%
%   网格布局:
%   - 交错网格: U速度在垂直面上(Ny+2 × Nx+1)
%   - 包括上下边界条件的虚拟单元
%   - 入口/出口边界通过狄利克雷/诺伊曼条件处理
%
%   See also INIT_GRID_V, INIT_GRID_P, GET_AB_U.

global Nx Ny h H L slip_opt

%% U速度的交错网格几何
% U速度位于压力单元的垂直面上
Ny_u = Ny + 2;                       % 包括上下虚拟行
Nx_u = Nx + 1;                       % 每行的垂直面数

% U速度节点的物理坐标
xx = h * (0:Nx);                     % 面中心: x = 0, h, 2h, ...
yy = h * (-0.5 : 1 : Ny + 0.5);      % 包括虚拟单元: y = -h/2, h/2, ...
[XX, YY] = meshgrid(xx, yy);

%% 基于几何的节点分类
% 从障碍物几何的主要固体/流体分类
is_solid = get_blocked_mask(block_info, XX, YY);
is_fluid = ~is_solid;

% 浸没边界检测(与固体相邻的流体节点)
is_solid_s = false(Ny_u, Nx_u); is_solid_s(2:end, :) = is_solid(1:end-1, :);
is_solid_n = false(Ny_u, Nx_u); is_solid_n(1:end-1, :) = is_solid(2:end, :);
is_solid_w = false(Ny_u, Nx_u); is_solid_w(:, 2:end) = is_solid(:, 1:end-1);
is_solid_e = false(Ny_u, Nx_u); is_solid_e(:, 1:end-1) = is_solid(:, 2:end);
is_solid_boundary = is_fluid & (is_solid_s | is_solid_n | is_solid_w | is_solid_e);

% 计算域边界识别
is_inlet        = false(Ny_u, Nx_u); is_inlet(2:end-1, 1) = true;      % 左边界
is_outlet       = false(Ny_u, Nx_u); is_outlet(2:end-1, end) = true;   % 右边界
is_ghost_bottom = false(Ny_u, Nx_u); is_ghost_bottom(1, :) = true;     % 底部虚拟
is_ghost_top    = false(Ny_u, Nx_u); is_ghost_top(end, :) = true;      % 顶部虚拟

%% 每个节点的方程类型分配
% 狄利克雷节点: 规定的速度值
is_dirichlet = is_solid | is_solid_boundary | is_inlet;

% 虚拟节点: 壁面边界条件(滑移/无滑移)
is_ghost = is_ghost_bottom | is_ghost_top;

% 诺伊曼节点: 零梯度出流条件
is_neumann = is_outlet;

% PDE节点: 求解完整动量方程
is_pde_solve = is_fluid & ~is_dirichlet & ~is_ghost & ~is_neumann;

%% 离散格式选择
% 边界附近使用1阶，内部使用2阶
is_not_pde_solve = ~is_pde_solve;
is_neighbor_s = false(Ny_u, Nx_u); is_neighbor_s(2:end, :) = is_not_pde_solve(1:end-1, :);
is_neighbor_n = false(Ny_u, Nx_u); is_neighbor_n(1:end-1, :) = is_not_pde_solve(2:end, :);
is_neighbor_w = false(Ny_u, Nx_u); is_neighbor_w(:, 2:end) = is_not_pde_solve(:, 1:end-1);
is_neighbor_e = false(Ny_u, Nx_u); is_neighbor_e(:, 1:end-1) = is_not_pde_solve(:, 2:end);
is_near_boundary = is_neighbor_s | is_neighbor_n | is_neighbor_w | is_neighbor_e;

is_1st_order = is_pde_solve & is_near_boundary;   % 边界附近保守
is_2nd_order = is_pde_solve & ~is_near_boundary;  % 内部精确

%% 常数矩阵系数的预组装
% 为稀疏矩阵组装创建全局节点索引
[MM, NN] = meshgrid(1:Nx_u, 1:Ny_u);
p_idx = (NN-1)*Nx_u + MM;                        % 当前节点索引
s_idx = circshift(p_idx, [1, 0]);                % 南邻居
n_idx = circshift(p_idx, [-1, 0]);               % 北邻居
w_idx = circshift(p_idx, [0, 1]);                % 西邻居

% 为稀疏矩阵组装预分配三元组数组
nnz_estimate = nnz(is_dirichlet) + 2*nnz(is_ghost) + 2*nnz(is_neumann);
I_const = zeros(nnz_estimate, 1);
J_const = zeros(nnz_estimate, 1);
V_const = zeros(nnz_estimate, 1);
current_pos = 0;

% 狄利克雷边界条件: u = 规定值
[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_dirichlet, ...
    {p_idx}, 1, ...
    I_const, J_const, V_const, current_pos, p_idx);

% 壁面边界条件(虚拟单元)
if strcmp(slip_opt, '无滑移')
    % 无滑移: (u_ghost + u_interior)/2 = 0
    [I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_bottom, ...
        {p_idx, n_idx}, [0.5; 0.5], ...
        I_const, J_const, V_const, current_pos, p_idx);
    [I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_top, ...
        {p_idx, s_idx}, [0.5; 0.5], ...
        I_const, J_const, V_const, current_pos, p_idx);
elseif strcmp(slip_opt, '滑移')
    % 滑移: u_ghost = u_interior
    [I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_bottom, ...
        {p_idx, n_idx}, [1; -1], ...
        I_const, J_const, V_const, current_pos, p_idx);
    [I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_ghost_top, ...
        {p_idx, s_idx}, [1; -1], ...
        I_const, J_const, V_const, current_pos, p_idx);
else
    error('无效的壁面条件选项。请选择"滑移"或"无滑移"。');
end

% 诺伊曼出流: ∂u/∂x = 0 → u_outlet = u_neighbor
[I_const, J_const, V_const, current_pos] = assemble_stencil_const(is_outlet, ...
    {p_idx, w_idx}, [1; -1], ...
    I_const, J_const, V_const, current_pos, p_idx);

% 移除未使用的预分配条目
I_const = I_const(1:current_pos);
J_const = J_const(1:current_pos);
V_const = V_const(1:current_pos);

%% 打包完整的网格结构
grid_u.h = h;
grid_u.L = L;
grid_u.H = H;
grid_u.Nx = Nx_u;
grid_u.Ny = Ny_u;
grid_u.xx = xx;
grid_u.yy = yy;

% 节点分类掩码
grid_u.is_solid = is_solid;
grid_u.is_solid_boundary = is_solid_boundary;
grid_u.is_inlet = is_inlet;
grid_u.is_outlet = is_outlet;
grid_u.is_ghost_bottom = is_ghost_bottom;
grid_u.is_ghost_top = is_ghost_top;

% 离散格式分类
grid_u.is_1st_order = is_1st_order;
grid_u.is_2nd_order = is_2nd_order;

% 预组装的常数矩阵系数
grid_u.I_const = I_const;
grid_u.J_const = J_const;
grid_u.V_const = V_const;

% 用于模板组装的节点索引映射
grid_u.p_idx = p_idx;
grid_u.s_idx = s_idx;
grid_u.n_idx = n_idx;
grid_u.w_idx = w_idx;
grid_u.e_idx = circshift(p_idx, [0, -1]);
grid_u.ss_idx = circshift(p_idx, [2, 0]);
grid_u.nn_idx = circshift(p_idx, [-2, 0]);
grid_u.ww_idx = circshift(p_idx, [0, 2]);
grid_u.ee_idx = circshift(p_idx, [0, -2]);

end
