function is_blocked = get_blocked_mask(block_info, XX, YY)
%GET_BLOCKED_MASK 为障碍物内部的网格点创建布尔掩码
%   此函数确定网格的哪些点位于一个或多个多边形障碍物内部。
%   这是设置复杂几何形状仿真的关键工具函数。
%
%   用法:
%       is_blocked = get_blocked_mask(block_info, XX, YY);
%
%   输入参数:
%       block_info - 元胞数组，每个元胞包含一个障碍物的信息。
%                    必须包含字段'.points'，这是一个M×2矩阵，
%                    包含定义障碍物多边形顶点的[x, y]坐标。
%       XX, YY     - 网格点的坐标矩阵，通常由MESHGRID生成。
%
%   输出参数:
%       is_blocked - 与XX和YY相同大小的布尔矩阵。如果网格点
%                    (XX(i,j), YY(i,j))位于任何指定障碍物的内部
%                    或边界上，则is_blocked(i,j)为true。否则为false。

% 获取网格的尺寸
[Ny, Nx] = size(XX);

% 将输出掩码初始化为全false
is_blocked = false(Ny, Nx);

% 将网格坐标矩阵展平为列向量，以便inpolygon函数高效处理
query_points_x = XX(:);
query_points_y = YY(:);

% 迭代处理block_info元胞数组中定义的每个障碍物
for k = 1:length(block_info)
    obstacle = block_info{k};
    poly_points = obstacle.points;
    
    % 使用MATLAB内置的inpolygon函数高效确定
    % 哪些查询点位于当前障碍物多边形的内部或边界上
    [in, ~] = inpolygon(query_points_x, query_points_y, poly_points(:,1), poly_points(:,2));
    
    % 将结果逻辑向量重新塑形为与原始网格尺寸匹配的2D掩码
    in_mask = reshape(in, Ny, Nx);
    
    % 使用逻辑OR运算将当前障碍物的掩码与主掩码合并
    is_blocked = is_blocked | in_mask;
end

end
