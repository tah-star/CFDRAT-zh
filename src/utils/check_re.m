% check_reynolds_number.m
function check_re(rho, U_char, mu, H, block_info, Re_limit)
%CHECK_RE 验证雷诺数在安全运行范围内
%   根据流动特性和计算域几何形状计算雷诺数。如果Re超过指定限制，
%   则终止仿真并显示详细错误信息。
%
%   输入参数:
%       rho        - 流体密度 [kg/m³]
%       U_char     - 特征速度 [m/s]
%       mu         - 动力粘度 [Pa·s]
%       H          - 计算域高度 [m]
%       block_info - 来自init_block的障碍物几何数据
%       Re_limit   - 最大允许雷诺数
%
%   特征长度选择规则:
%   - 无障碍物：使用计算域高度H
%   - 有障碍物：使用最大障碍物高度
%   - 确保给定离散化方案的数值稳定性
%
%   See also INIT_BLOCK, SETUP_SIMULATION.

%% 确定特征长度尺度
if isempty(block_info)
    % 纯通道流：特征长度为计算域高度
    L_char = H;
else
    % 含障碍物流动：使用最大障碍物高度
    max_obstacle_height = 0;
    for i = 1:numel(block_info)
        y_min_obstacle = block_info{i}.y_coords(1);
        y_max_obstacle = block_info{i}.y_coords(2);
        current_obstacle_height = y_max_obstacle - y_min_obstacle;
        if current_obstacle_height > max_obstacle_height
            max_obstacle_height = current_obstacle_height;
        end
    end
    L_char = max_obstacle_height;
end

% 处理退化情况（薄障碍物）
if L_char <= 1e-9
    L_char = H;  % 回退到计算域高度
end

%% 计算雷诺数
if mu > 0
    Re = rho * U_char * L_char / mu;
else
    Re = inf;  % 无粘流极限
end

%% 验证是否超过稳定性限制
if Re > Re_limit
    % 为超出雷诺数创建详细错误信息
    error_message = sprintf([...
        '雷诺数超出安全运行限制。\n\n' ...
        '  - 计算得到的Re: %.0f\n' ...
        '  - 最大允许值: %.0f\n\n' ...
        '为保持数值稳定性，请考虑:\n' ...
        '  - 降低入口速度\n' ...
        '  - 降低流体密度\n' ...
        '  - 增加流体粘度'], ...
        Re, Re_limit);
    
    error('MyApp:ReynoldsNumberExceeded', error_message);
end

% 如果雷诺数可接受，函数完成

end
