function [I, J, V, current_pos] = assemble_stencil_const(mask, stencil_J_maps, stencil_V, I, J, V, current_pos, p_idx_map)
%ASSEMBLE_STENCIL_CONST 为稀疏矩阵组装填充IJV三元组（常系数模板）
%   将有限差分模板高效地组装成稀疏矩阵三元组格式。处理由mask标识的
%   节点，应用给定的模板模式，具有高性能特征。
%
%   输入参数:
%       mask           - 标识目标节点的逻辑矩阵
%       stencil_J_maps - 模板邻节点索引映射的元胞数组
%       stencil_V      - 模板系数值的列向量
%       I, J, V        - 正在填充的三元组数组
%       current_pos    - 三元组数组中的当前位置
%       p_idx_map      - 中心节点的网格到线性索引映射
%
%   输出参数:
%       I, J, V        - 更新后的三元组数组
%       current_pos    - 插入后的新位置
%
%   算法说明:
%   对于N个节点和M点模板，使用交错结构生成N×M个三元组，
%   确保正确的稀疏矩阵组装。
%
%   示例：5点模板创建如下模式：
%   I: [p1,p2,...,pN, p1,p2,...,pN, ...]  (重复M次)
%   J: [s1,s2,...,sN, w1,w2,...,wN, ...]  (邻节点索引)
%   V: [vs,vs,...,vs, vw,vw,...,vw, ...]  (系数值)
%
%   See also SPARSE, ASSEMBLE_STENCIL_DYNAMIC.

% 初始化并检查空mask
len = nnz(mask);
if len == 0
    return; 
end

num_stencil_pts = numel(stencil_J_maps);
num_new_entries = num_stencil_pts * len;
indices = current_pos+1 : current_pos+num_new_entries;

% 提取中心点的线性索引
p_indices_masked = p_idx_map(mask);

% 填充行索引(I): 每个模板点重复中心点索引
I(indices) = repmat(p_indices_masked, num_stencil_pts, 1);

% 填充列索引(J): 按模板位置组织的邻节点索引
J_block = zeros(num_new_entries, 1);
for k = 1:num_stencil_pts
    neighbor_map = stencil_J_maps{k};
    J_k_indices = neighbor_map(mask);
    
    start_idx = (k-1)*len + 1;
    end_idx = k*len;
    J_block(start_idx:end_idx) = J_k_indices;
end
J(indices) = J_block;

% 填充系数值(V): 每个节点重复模板值
V_block = zeros(num_new_entries, 1);
for k = 1:num_stencil_pts
    value_k = stencil_V(k);
    
    start_idx = (k-1)*len + 1;
    end_idx = k*len;
    V_block(start_idx:end_idx) = value_k;
end
V(indices) = V_block;
current_pos = current_pos + num_new_entries;
end
