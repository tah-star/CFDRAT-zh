function run_simu(varargin)
%RUN_SIMU 2D不可压缩流动仿真的主驱动程序
%   RUN_SIMU(PARAMS)在交错网格上使用PISO算法执行完整的CFD仿真。
%   处理通过二值图像定义的复杂几何形状。
%
%   RUN_SIMU() 使用内置默认参数运行，用于快速测试。
%   RUN_SIMU(PARAMS) 使用来自GUI或脚本的用户指定参数运行。
%
%   输入参数:
%       PARAMS - (可选)包含仿真参数的结构:
%                .H                 - 计算域高度 [m]
%                .u_inlet_func_str  - 入口速度表达式字符串
%                .t_simu           - 总仿真时间 [s]
%                .T_record         - 数据输出间隔 [s]
%                .rho              - 流体密度 [kg/m³]
%                .mu               - 动力粘度 [Pa*s]
%                .node_scale       - 近似总网格节点数
%                .speed_opt        - 基于CFL的时间步长 ('fast'|'medium')
%                .slip_opt         - 壁面条件 ('no-slip'|'slip')
%                .ramp_opt         - 入口渐变类型 ('linear'|'smoothstep'|'none')
%                .image_filename   - 障碍物图像文件路径
%                .save_filename    - 输出.mat文件路径
%                .image_display_opt- 可视化模式
%
%   示例:
%       % 使用默认值快速测试
%       run_simu();
%       
%       % 自定义仿真
%       params.H = 0.5;
%       params.t_simu = 2.0;
%       params.image_filename = 'cylinder.png';
%       run_simu(params);
%
%   仿真流水线:
%   1. 初始化几何形状和网格
%   2. 检查雷诺数稳定性
%   3. 设置线性求解器和预条件器
%   4. 执行时间步进循环(PISO算法)
%   5. 保存结果并启动可视化
%
%   See also CFDRAT, SOLVE_PDE, INIT_BLOCK, DRAW_UV.

% 解析输入参数 - 如果未提供则使用默认值
if nargin == 0 || isempty(varargin{1})
    % 默认测试案例: 中等雷诺数下绕障碍物流动
    params = struct();
    params.H = 0.3;                                 % 计算域高度
    params.u_inlet_func_str = '0.15*4*y*(H-y)/H^2'; % 抛物线入口剖面
    params.t_ramp = 0.5;                            % 入口渐变时间
    params.slip_opt = "无滑移";                    % 壁面边界条件
    params.ramp_opt = "线性上升";                     % 入口渐变类型
    params.rho = 1;                                 % 流体密度
    params.mu = 1e-4;                              % 动力粘度
    params.t_simu = 1;                             % 总仿真时间
    params.T_record = 0.05;                        % 输出时间间隔
    params.node_scale = 5e4;                       % 目标网格分辨率
    params.speed_opt = "快";                     % 基于CFL的时间步进
    params.image_filename = "test_fig.png";        % 障碍物几何文件
    params.save_filename = "sim_result.mat";       % 结果输出文件
    params.image_display_opt = "内部色块";       % 可视化偏好
    params.enable_parallel = '关';
else
    % 使用来自GUI或外部脚本的参数
    params = varargin{1};
end

% 初始化全局仿真参数
global H L Nx Ny h dt rho mu slip_opt simulation_control

enable_parallel = params.enable_parallel;               % 并行计算标志
mu = params.mu;                                   % 提取流体属性
rho = params.rho;
t_simu = params.t_simu;                          % 提取时间参数
T_record = params.T_record;
node_scale = params.node_scale;                  % 提取网格参数
H = params.H;                                    % 计算域高度
u_inlet_func_str = params.u_inlet_func_str;     % 提取边界条件
t_ramp = params.t_ramp;                          % 入口速度渐变时间
slip_opt = params.slip_opt;                      % 壁面边界条件类型
ramp_opt = params.ramp_opt;                      % 入口渐变函数类型
speed_opt = params.speed_opt;                    % 提取求解器选项
image_filename = params.image_filename;          % 提取I/O参数
image_display_opt = params.image_display_opt;   % 后处理可视化模式
save_filename = params.save_filename;            % 结果输出文件路径

% 数值求解器配置(为稳定性精细调节)
Re_limit = 300;                                  % 最大雷诺数
threshold_increasment_ratio = 1.6;               % 预条件器更新触发器
threshold_maxmin_ratio = 2.5;                    % 最小/最大迭代比率触发器
uplimit_recorder = 20;                           % 预条件器重用计数上限
max_iter = 200;                                  % 线性求解器迭代限制
residual_tol = 5e-4;                            % 收敛容差

% 执行完整的仿真流水线
u_inlet_func = get_inlet_func(u_inlet_func_str, ramp_opt, t_ramp, H, slip_opt);  % 创建入口速度函数句柄
[u_inlet_max, u_all_max] = get_u_max(u_inlet_func, H, t_simu, t_ramp);         % 为CFL估算最大速度
[block_info, L, h, Nx, Ny] = init_block(image_filename, H, node_scale);        % 生成几何和网格参数
check_re(rho, u_inlet_max, mu, H, block_info, Re_limit);                       % 验证雷诺数稳定性
dt = init_dt(speed_opt, u_all_max);                                            % 计算稳定时间步长
grid_u = init_grid_u(block_info);                                              % 初始化速度网格
grid_v = init_grid_v(block_info);                                              % 初始化v速度网格和边界条件矩阵
[A_dp, grid_p] = init_grid_p(block_info);                                      % 初始化压力网格和矩阵
continue_flag = confirm_grid(grid_p, block_info);                              % 可视化网格验证
if ~continue_flag; return; end
precond = init_precond(threshold_increasment_ratio, threshold_maxmin_ratio, uplimit_recorder, A_dp);  % 设置自适应预条件器
solver_opt = init_solver(max_iter, residual_tol);                              % 配置迭代求解器参数
init_parpool(enable_parallel);                                                 % 设置并行计算
[u_all, v_all, p_all, t_solve] = solve_pde(t_simu, T_record, u_inlet_func, grid_u, grid_v, grid_p, A_dp, precond, solver_opt);  % 执行PISO时间步进循环
if strcmp(simulation_control, 'stopped')
    return;
else
    save_data(u_all, v_all, p_all, u_inlet_func, node_scale, dt, t_simu, T_record, t_solve, grid_u, grid_v, grid_p, block_info, save_filename);  % 打包和保存结果
    draw_uv(save_filename, image_display_opt);                                     % 启动交互式可视化
end

end
