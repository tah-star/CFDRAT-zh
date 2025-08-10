function [A_dp, grid_p] = init_grid_p(block_info)
%INIT_GRID_P 初始化压力网格并组装泊松矩阵
%   创建压力场网格结构并组装压力泊松方程∇²p = RHS的常系数矩阵。
%   处理具有全面边界条件分类的复杂几何形状。
%
%   输入参数:
%       block_info - 障碍物边界信息的元胞数组
%
%   输出参数:
%       A_dp    - 稀疏压力泊松矩阵(正定)
%       grid_p  - 具有节点分类的完整压力网格结构
%
%   主要特性:
%   - 互斥的节点分类(内部、边界、角点等)
%   - 带边界修正的5点有限差分模板
%   - 稀疏矩阵组装的高效IJV三元组格式
%   - 处理诺伊曼(壁面)和狄利克雷(出口)边界条件
%
%   节点类型:
%   - 内部: 标准5点拉普拉斯算子
%   - 单一边界: 4点模板
%   - 角点: 3点模板
%   - 凹角: 2点模板
%
%   See also INIT_GRID_U, INIT_GRID_V, CONFIRM_GRID.

global Nx Ny h H L

%% 压力网格几何(单元中心)
Ny_p = Ny;
Nx_p = Nx;
NyNx_p = Ny_p * Nx_p;                            % 压力节点总数

% 单元中心的物理坐标
xx = h * (0.5 : 1 : Nx-0.5);                    % x = h/2, 3h/2, ...
yy = h * (0.5 : 1 : Ny-0.5);                    % y = h/2, 3h/2, ...
[XX, YY] = meshgrid(xx, yy);

%% 主要几何分类
% 从障碍物几何确定固体/流体状态
is_solid = get_blocked_mask(block_info, XX, YY);
is_fluid = ~is_solid;

% 检测浸没边界(与固体相邻的流体节点)
is_solid_s = false(Ny_p, Nx_p); is_solid_s(2:end, :) = is_solid(1:end-1, :); % 南邻是固体
is_solid_n = false(Ny_p, Nx_p); is_solid_n(1:end-1, :) = is_solid(2:end, :); % 北邻是固体
is_solid_w = false(Ny_p, Nx_p); is_solid_w(:, 2:end) = is_solid(:, 1:end-1); % 西邻是固体
is_solid_e = false(Ny_p, Nx_p); is_solid_e(:, 1:end-1) = is_solid(:, 2:end); % 东邻是固体
is_solid_boundary = is_fluid & (is_solid_s | is_solid_w | is_solid_e | is_solid_n);

% 计算域边界
is_domain_b = false(Ny_p, Nx_p); is_domain_b(1, :) = true;     % 底壁
is_domain_t = false(Ny_p, Nx_p); is_domain_t(end, :) = true;   % 顶壁
is_domain_l = false(Ny_p, Nx_p); is_domain_l(:, 1) = true;     % 左壁(入口)
is_domain_r = false(Ny_p, Nx_p); is_domain_r(:, end) = true;   % 右壁(出口)

%% 边界条件类型识别
% 诺伊曼边界: ∂p/∂n = 0 (壁面和浸没边界)
is_neumann_s = (is_domain_b | (is_fluid & is_solid_s)) & ~is_solid;
is_neumann_n = (is_domain_t | (is_fluid & is_solid_n)) & ~is_solid;
is_neumann_w = (is_domain_l | (is_fluid & is_solid_w)) & ~is_solid;
is_neumann_e = (is_fluid & is_solid_e) & ~is_solid;
% 狄利克雷边界: p = 0 (出口参考压力)
is_dirichlet_e = is_domain_r & ~is_neumann_e & ~is_solid;

%% 互斥节点分类
% 计算每个节点的边界条件数以进行系统分类
num_neumann_boundaries = double(is_neumann_s) + double(is_neumann_n) + double(is_neumann_w) + double(is_neumann_e);
num_dirichlet_boundaries = double(is_dirichlet_e);

% 内部流体节点(标准5点模板)
is_fluid_interior = is_fluid & (num_neumann_boundaries == 0) & (num_dirichlet_boundaries == 0);

% 单一边界节点(4点模板)
is_neumann_s_only = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 0) & is_neumann_s;
is_neumann_n_only = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 0) & is_neumann_n;
is_neumann_w_only = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 0) & is_neumann_w;
is_neumann_e_only = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 0) & is_neumann_e;
is_dirichlet_e_only = is_fluid & (num_dirichlet_boundaries == 1) & (num_neumann_boundaries == 0) & is_dirichlet_e;

% 具有两个边界的角点节点(3点模板)
is_corner_S_neumann_W_neumann = is_fluid & (num_neumann_boundaries == 2) & (num_dirichlet_boundaries == 0) & is_neumann_s & is_neumann_w;
is_corner_N_neumann_W_neumann = is_fluid & (num_neumann_boundaries == 2) & (num_dirichlet_boundaries == 0) & is_neumann_n & is_neumann_w;
is_corner_S_neumann_E_neumann = is_fluid & (num_neumann_boundaries == 2) & (num_dirichlet_boundaries == 0) & is_neumann_s & is_neumann_e;
is_corner_N_neumann_E_neumann = is_fluid & (num_neumann_boundaries == 2) & (num_dirichlet_boundaries == 0) & is_neumann_n & is_neumann_e;
% 混合角点(诺伊曼+狄利克雷)
is_corner_S_neumann_E_dirichlet = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 1) & is_neumann_s & is_dirichlet_e;
is_corner_N_neumann_E_dirichlet = is_fluid & (num_neumann_boundaries == 1) & (num_dirichlet_boundaries == 1) & is_neumann_n & is_dirichlet_e;

% 具有三个边界的凹角节点(2点模板)
is_concave_NWS = is_fluid & (num_neumann_boundaries == 3) & (num_dirichlet_boundaries == 0) & is_neumann_n & is_neumann_w & is_neumann_s;
is_concave_NES = is_fluid & (num_neumann_boundaries == 3) & (num_dirichlet_boundaries == 0) & is_neumann_n & is_neumann_e & is_neumann_s;
is_concave_WNE = is_fluid & (num_neumann_boundaries == 3) & (num_dirichlet_boundaries == 0) & is_neumann_w & is_neumann_n & is_neumann_e;
is_concave_WSE = is_fluid & (num_neumann_boundaries == 3) & (num_dirichlet_boundaries == 0) & is_neumann_w & is_neumann_s & is_neumann_e;

%% 使用IJV三元组格式的稀疏矩阵组装
% 为高效模板组装创建节点索引映射
[II, JJ] = meshgrid(1:Nx_p, 1:Ny_p);
p = (JJ-1)*Nx_p + II;                            % 当前节点索引
s = circshift(p, [1, 0]); n = circshift(p, [-1, 0]); % 南、北邻居
w = circshift(p, [0, 1]); e = circshift(p, [0, -1]); % 西、东邻居

% 预分配三元组数组
nnz_estimate = nnz(is_solid) + 5*nnz(is_fluid);  % 保守估计
I = zeros(nnz_estimate, 1); J = zeros(nnz_estimate, 1); V = zeros(nnz_estimate, 1);
current_pos = 0;

% 为每种节点分类组装模板系数
% 固体节点: 恒等方程(确保非奇异性)
[I, J, V, current_pos] = assemble_stencil_const(is_solid, {p}, -1, I, J, V, current_pos, p);
% 内部: 标准5点拉普拉斯算子
[I, J, V, current_pos] = assemble_stencil_const(is_fluid_interior, {s,w,p,e,n}, [1;1;-4;1;1], I, J, V, current_pos, p);
% 单一边界节点: 修正的4点模板
[I, J, V, current_pos] = assemble_stencil_const(is_neumann_s_only, {w,p,e,n}, [1;-3;1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_neumann_n_only, {s,w,p,e}, [1;1;-3;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_neumann_w_only, {s,p,e,n}, [1;-3;1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_neumann_e_only, {s,w,p,n}, [1;1;-3;1], I, J, V, current_pos, p);
% 角点节点: 3点模板
[I, J, V, current_pos] = assemble_stencil_const(is_corner_S_neumann_W_neumann, {p,e,n}, [-2;1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_corner_N_neumann_W_neumann, {s,p,e}, [1;-2;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_corner_S_neumann_E_neumann, {w,p,n}, [1;-2;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_corner_N_neumann_E_neumann, {s,w,p}, [1;1;-2], I, J, V, current_pos, p);
% 凹角: 2点模板
[I, J, V, current_pos] = assemble_stencil_const(is_concave_NWS, {p,e}, [-1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_concave_NES, {p,w}, [-1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_concave_WNE, {p,s}, [-1;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_concave_WSE, {p,n}, [-1;1], I, J, V, current_pos, p);
% 狄利克雷边界: 参考压力的增强模板
[I, J, V, current_pos] = assemble_stencil_const(is_dirichlet_e_only, {s,w,p,n}, [1;1;-5;1], I, J, V, current_pos, p);
% 混合边界角点
[I, J, V, current_pos] = assemble_stencil_const(is_corner_S_neumann_E_dirichlet, {w,p,n}, [1;-4;1], I, J, V, current_pos, p);
[I, J, V, current_pos] = assemble_stencil_const(is_corner_N_neumann_E_dirichlet, {s,w,p}, [1;1;-4], I, J, V, current_pos, p);

% 完成稀疏矩阵
I = I(1:current_pos); J = J(1:current_pos); V = V(1:current_pos);
V = -V / h^2;                                    % 符号翻转以获得正定性 + 网格缩放
A_dp = sparse(I, J, V, NyNx_p, NyNx_p);

%% 打包完整的网格结构
grid_p.h = h; grid_p.L = L; grid_p.H = H;
grid_p.Nx_p = Nx_p; grid_p.Ny_p = Ny_p;
grid_p.xx = xx; grid_p.yy = yy;
grid_p.is_solid = is_solid; grid_p.is_fluid = is_fluid;
grid_p.is_solid_boundary = is_solid_boundary;
grid_p.is_fluid_interior = is_fluid_interior;
grid_p.is_neumann_s_only = is_neumann_s_only;
grid_p.is_neumann_n_only = is_neumann_n_only;
grid_p.is_neumann_w_only = is_neumann_w_only;
grid_p.is_neumann_e_only = is_neumann_e_only;
grid_p.is_dirichlet_e_only = is_dirichlet_e_only;
grid_p.is_corner_S_neumann_W_neumann = is_corner_S_neumann_W_neumann;
grid_p.is_corner_N_neumann_W_neumann = is_corner_N_neumann_W_neumann;
grid_p.is_corner_S_neumann_E_neumann = is_corner_S_neumann_E_neumann;
grid_p.is_corner_N_neumann_E_neumann = is_corner_N_neumann_E_neumann;
grid_p.is_corner_S_neumann_E_dirichlet = is_corner_S_neumann_E_dirichlet;
grid_p.is_corner_N_neumann_E_dirichlet = is_corner_N_neumann_E_dirichlet;
grid_p.is_concave_NWS = is_concave_NWS;
grid_p.is_concave_NES = is_concave_NES;
grid_p.is_concave_WNE = is_concave_WNE;
grid_p.is_concave_WSE = is_concave_WSE;

end
