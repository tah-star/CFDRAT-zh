function save_filename = save_data(u_all, v_all, p_all, u_inlet_func, node_scale, dt, t_simu, T_record, t_solve, grid_u, grid_v, grid_p, block_info, save_filename)
%SAVE_DATA 将CFD仿真结果处理并保存到MAT文件
%   通过将交错网格速度场插值到一致的网格位置来后处理仿真数据，
%
%   输入参数:
%       u_all, v_all, p_all - 速度和压力场的时间序列
%       u_inlet_func        - 入口边界条件函数句柄
%       node_scale, dt      - 网格分辨率和时间步参数
%       t_simu, T_record    - 仿真和记录时间间隔
%       t_solve            - 总计算时间 [s]
%       grid_u/v/p         - 交错网格的网格结构
%       block_info         - 计算域几何信息
%       save_filename      - 输出文件路径
%
%   输出参数:
%       save_filename - 确认的输出文件路径
%
%   处理过程:
%   - 将U/V速度从面心插值到一致的位置
%   - 将压力从单元中心插值到网格节点
%   - 转换为单精度以提高存储效率
%   - 打包仿真参数和元数据
%   - 以压缩MAT格式保存以兼容大数据集
%
%   另见 LOAD, INTERP_STAGGERED_TO_NODES.

%% 确保输出目录存在
[save_path, ~, ~] = fileparts(save_filename);
if ~isempty(save_path) && ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('Created output directory: %s\n', save_filename);
end

%% 初始化数据结构并分配内存
info = struct();

[Ny_u, Nx_u, num_frames] = size(u_all);
[Ny_v, Nx_v, ~] = size(v_all);
[Ny_p, Nx_p, ~] = size(p_all);

% 用最优数据类型预分配存储空间
info.u_all = zeros(Ny_u - 1, Nx_u, num_frames, 'single');
info.v_all = zeros(Ny_v, Nx_v - 1, num_frames, 'single');
info.p_all = zeros(Ny_p + 1, Nx_p + 1, num_frames, 'single');  % 节点处的压力

%% 通过插值处理速度和压力场
for k = 1:num_frames
    % U速度: 从垂直面心插值到角点节点
    u_slice = single(u_all(:,:,k));
    info.u_all(:,:,k) = (u_slice(1:end-1,:) + u_slice(2:end,:)) / 2;
    
    % V速度: 从水平面心插值到角点节点
    v_slice = single(v_all(:,:,k));
    info.v_all(:,:,k) = (v_slice(:,1:end-1) + v_slice(:,2:end)) / 2;
    
    % 压力: 从单元中心插值到网格节点
    p_slice = single(p_all(:,:,k));
    p_nodes = interp_p(p_slice);
    info.p_all(:,:,k) = p_nodes;
end

h = grid_u.h;
xx = h*(0:Nx_p);
yy = h*(0:Ny_p);
[XX, YY] = meshgrid(xx, yy);

%% 打包仿真元数据和参数
info.u_inlet_func = u_inlet_func;
info.dt = dt;
info.grid_u = grid_u;
info.grid_v = grid_v;
info.grid_p = grid_p;
info.h = h;
info.block_info = block_info;
info.T_record = T_record;
info.t_simu = t_simu;
info.node_scale = node_scale;
info.t_solve = t_solve;
info.XX = XX;
info.YY = YY;
info.save_time = datetime('now');
 
%% 将数据写入压缩MAT文件
% 使用v7.3格式支持大文件和压缩
save(save_filename, 'info', '-v7.3');
fprintf('Data successfully saved to: %s\n', save_filename);

end


function p_nodes_all = interp_p(p_centers_all)
%   使用向量化操作将压力场从(Ny x Nx)单元中心转换为
%   (Ny+1 x Nx+1)网格节点
%   边界条件:
%   - 上/左/下边界: 零梯度 (∂p/∂n = 0)  
%   - 右边界: 零值 (p = 0)
%
%   输入:
%       p_centers_all - 单元中心处的压力 [Ny x Nx]
%   输出:
%       p_nodes_all   - 网格节点处的压力 [Ny+1 x Nx+1]

[Ny, Nx] = size(p_centers_all);
p_nodes_all = zeros(Ny + 1, Nx + 1, 'single');

% 内部节点 (2:Ny, 2:Nx) 使用周围4个单元中心的平均值
p_nodes_all(2:Ny, 2:Nx) = 0.25 * (...
    p_centers_all(1:Ny-1, 1:Nx-1) + ...  % 左上
    p_centers_all(1:Ny-1, 2:Nx) + ...  % 右上
    p_centers_all(2:Ny,   1:Nx-1) + ...  % 左下
    p_centers_all(2:Ny,   2:Nx));      % 右下

% 下边界 (i=1): 零梯度 ∂p/∂y = 0
p_nodes_all(1, 2:Nx) = p_nodes_all(2, 2:Nx);

% 上边界 (i=Ny+1): 零梯度 ∂p/∂y = 0
p_nodes_all(Ny+1, 2:Nx) = p_nodes_all(Ny, 2:Nx);

% 左边界 (j=1): 零梯度 ∂p/∂x = 0  
p_nodes_all(2:Ny, 1) = p_nodes_all(2:Ny, 2);

% 右边界 (j=Nx+1): 零值 p = 0
p_nodes_all(:, Nx+1) = 0;

% 左下角: 两个方向的零梯度
p_nodes_all(1, 1) = p_nodes_all(2, 2);

% 左上角: 两个方向的零梯度
p_nodes_all(Ny+1, 1) = p_nodes_all(Ny, 2);

% 右下角和右上角: 零值 (右边界条件主导) 无需分配

end
