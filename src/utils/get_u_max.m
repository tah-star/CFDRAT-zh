% estimate_max_inlet_velocity.m
function [u_inlet_max, u_all_max] = get_u_max(inlet_func, H, t_end, t_ramp)
%GET_U_MAX 估算最大入口速度用于CFL和稳定性分析
%   在时间和空间上鲁棒地采样入口速度剖面以找到绝对最大速度。
%   对确定稳定时间步长和雷诺数计算至关重要。
%
%   输入参数:
%       inlet_func - 入口速度剖面的函数句柄@(t,y) [m/s]
%       H          - 计算域高度 [m]
%       t_end      - 总仿真时间 [s]
%       t_ramp     - 渐变持续时间 [s]
%
%   输出参数:
%       u_inlet_max - 最大入口速度幅值 [m/s]
%       u_all_max   - 域内最大速度的保守估计 [m/s]
%                     (假设由于加速效应有2倍放大)
%
%   方法:
%   - 在精细的时空网格上采样速度
%   - 聚焦采样到可能出现峰值的渐变完成附近
%   - 返回CFL计算的保守域内估计
%
%   See also INIT_DT, CHECK_RE.

% 定义采样分辨率以进行鲁棒的最大值检测
num_y_points = 101;                             
y_samples = linspace(0, H, num_y_points);       

% 聚焦时间采样到关键的渐变完成期
num_t_points = 51;
% 在可能出现最大值的渐变完成附近集中采样
t_samples = unique([0, linspace(t_ramp * 0.8, t_ramp * 1.2, num_t_points), t_end]);

% 为矢量化速度评估创建网格
[T_grid, Y_grid] = meshgrid(t_samples, y_samples);

% 在整个时空域上评估入口速度函数
velocity_values = inlet_func(T_grid, Y_grid);

% 从所有采样点提取绝对最大速度
u_inlet_max = max(abs(velocity_values(:)));

% 数值问题的安全检查
if isempty(u_inlet_max) || isnan(u_inlet_max) || isinf(u_inlet_max)
    u_inlet_max = 0;                          
end

u_all_max = 2 * u_inlet_max;                     % 域内最大值的保守估计

end
