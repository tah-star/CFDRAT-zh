% file: fill_IJV_dynamic.m
function [I, J, V, current_pos] = assemble_stencil_dynamic(mask, stencil_J_maps, stencil_C_maps, I, J, V, current_pos, p_idx_map)
%ASSEMBLE_STENCIL_DYNAMIC 为动态系数模板填充IJV三元组
%   类似于ASSEMBLE_STENCIL_CONST，但处理空间和时间变化的模板系数。
%   用于对流扩散算子，其系数取决于局部流动条件。
%
%   输入参数:
%       mask            - 标识目标节点的逻辑矩阵
%       stencil_J_maps  - 邻节点索引映射的元胞数组
%       stencil_C_maps  - 系数矩阵的元胞数组（在空间/时间中变化）
%       I, J, V         - 正在填充的三元组数组
%       current_pos     - 三元组数组中的当前位置
%       p_idx_map       - 中心节点的网格到线性索引映射
%
%   输出参数:
%       I, J, V         - 更新后的三元组数组
%       current_pos     - 插入后的新位置
%
%   与CONST版本的关键区别:
%   - stencil_C_maps包含矩阵，而不是标量
%   - 从这些矩阵中逐节点提取系数值
%
%   用途：动量方程组装，其中对流系数取决于局部速度场
%   并在每个时间步长中发生变化。
%
%   See also ASSEMBLE_STENCIL_CONST, GET_AB_U, GET_AB_V.

len = nnz(mask);
if len == 0, return; end

num_stencil_pts = numel(stencil_J_maps);
num_new_entries = num_stencil_pts * len;
indices = current_pos+1 : current_pos+num_new_entries;

% 行索引：每个模板位置重复中心点
I(indices) = repmat(p_idx_map(mask), num_stencil_pts, 1);

J_block = zeros(num_new_entries, 1);
V_block = zeros(num_new_entries, 1);
for k = 1:num_stencil_pts
    idx_range = (k-1)*len+1 : k*len;
    
    % 从邻节点映射获取列索引
    J_block(idx_range) = stencil_J_maps{k}(mask);
    
    % 从矩阵获取系数值（空间变化）
    V_block(idx_range) = stencil_C_maps{k}(mask);
end
J(indices) = J_block;
V(indices) = V_block;

current_pos = current_pos + num_new_entries;
end
