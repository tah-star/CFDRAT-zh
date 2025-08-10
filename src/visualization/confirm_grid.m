function continue_flag = confirm_grid(grid_p, block_info)
%CONFIRM_GRID 交互式网格可视化和完整性验证
%   显示带颜色编码的压力网格节点分类，并在仿真前进行全面完整性检查。
%   提供用户界面用于网格确认或参数调整。
%
%   输入参数:
%       grid_p     - 包含节点分类掩码的压力网格结构
%       block_info - 障碍物边界信息的元胞数组
%
%   输出参数:
%       continue_flag - 用户决定：true（继续）或false（返回设置）
%
%   功能特性:
%   - 所有节点类型的颜色编码可视化
%   - 覆盖检查：确保所有节点都被分类
%   - 唯一性检查：验证没有重叠分类
%   - 交互式缩放和平移以进行详细检查
%   - 用户确认界面
%
%   网格完整性测试:
%   - 覆盖性：每个节点必须属于且仅属于一个类别
%   - 唯一性：没有节点可以有多个分类
%   - 包括所有边界条件类型（角点、凹角等）
%
%   See also INIT_GRID_P, DRAW_GRID.

%% 初始化可视化窗口
fig_grid = figure('Name', '网格节点分类 - 检查并继续', ...
                  'NumberTitle', 'off', 'Position', [150, 150, 1400, 800]);
ax = gca;
hold(ax, 'on');
axis(ax, 'equal');
box(ax, 'on');

% 提取网格参数
Ny = grid_p.Ny_p; Nx = grid_p.Nx_p;
h = grid_p.h; L = grid_p.L; H = grid_p.H;
[XX, YY] = meshgrid(grid_p.xx, grid_p.yy);

fprintf('--- 网格概览 ---\n');
fprintf('计算域尺寸: %.3f m (L) x %.3f m (H)\n', L, H);
fprintf('网格分辨率: %d (Nx) x %d (Ny)\n', Nx, Ny);
fprintf('网格尺寸 (h): %.4f m\n', h);

%% 定义可视化类别以便于清晰显示
% 将计算域边界与其他分类分开
is_inlet = false(Ny, Nx); is_inlet(2:end-1, 1) = true;
is_inlet(grid_p.is_solid(:,1)) = false;            % 排除固体节点
is_outlet = false(Ny, Nx); is_outlet(2:end-1, end) = true;
is_outlet(grid_p.is_solid(:,end)) = false;         % 排除固体节点
is_walls = false(Ny, Nx); is_walls([1, end], :) = true;
is_walls(grid_p.is_solid([1, end],:)) = false;     % 排除固体节点
is_visual_fluid = grid_p.is_fluid & ~is_inlet & ~is_outlet & ~is_walls & ~grid_p.is_solid_boundary;

% 定义每个类别的绘制样式
plotting_categories = {
    is_visual_fluid,            '流体内部',                [0.2 0.8 0.6], '.', 8;
    is_inlet,                   '计算域入口',              [0.1 0.3 0.7], '.', 20;
    is_outlet,                  '计算域出口',              [0.7 0.1 0.3], '.', 20;
    is_walls,                   '计算域壁面(上/下)',        [0.2 0.2 0.2], '.', 20;
    grid_p.is_solid_boundary,   '浸没边界(固体)',          [0.9 0.4 0.1], '.', 20;
    grid_p.is_solid,            '固体内部(忽略)',          [0.8 0.8 0.8], 'x', 12;
};

%% 绘制带图例的网格分类
plot_handles = [];
legend_entries = {};

for i = 1:size(plotting_categories, 1)
    mask = plotting_categories{i, 1};
    label = plotting_categories{i, 2};
    color = plotting_categories{i, 3};
    marker = plotting_categories{i, 4};
    markersize = plotting_categories{i, 5};
    
    if any(mask(:))
        h_plot = plot(ax, XX(mask), YY(mask), marker, 'Color', color, ...
                      'MarkerSize', markersize, 'LineWidth', 1.3);
        plot_handles(end+1) = h_plot;
        legend_entries{end+1} = sprintf('%s (%d 节点)', label, nnz(mask));
    end
end

% 叠加障碍物边界作为参考
if ~isempty(block_info)
    boundary = block_info{1}.points;
    h_outline = plot(ax, [boundary(:,1); boundary(1,1)], [boundary(:,2); boundary(1,2)], 'k-', 'LineWidth', 2);
    plot_handles(end+1) = h_outline;
    legend_entries{end+1} = '真实障碍物轮廓';
    for k = 2:length(block_info)
        boundary = block_info{k}.points;
        plot(ax, [boundary(:,1); boundary(1,1)], [boundary(:,2); boundary(1,2)], 'k-', 'LineWidth', 2);
    end
end

%% 完善图形外观
axis(ax, 'tight');
title(ax, '网格预览: 可放大检查细节，完成后点击下方按钮', 'FontSize', 14, 'FontWeight', 'bold');
xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)');
set(ax, 'FontSize', 11, 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.5);
legend(plot_handles, legend_entries, 'Location', 'northeast', 'FontSize', 10, 'Box', 'on');
grid(ax, 'on');

%% 网格完整性验证
fprintf('\n--- 网格完整性检查 ---\n');

% 所有节点分类掩码的综合列表
all_masks_for_check = { 
    grid_p.is_solid; 
    grid_p.is_fluid_interior;
    grid_p.is_neumann_s_only; 
    grid_p.is_neumann_n_only; 
    grid_p.is_neumann_w_only; 
    grid_p.is_neumann_e_only;
    grid_p.is_dirichlet_e_only;
    grid_p.is_corner_S_neumann_W_neumann;
    grid_p.is_corner_N_neumann_W_neumann;
    grid_p.is_corner_S_neumann_E_neumann;
    grid_p.is_corner_N_neumann_E_neumann;
    grid_p.is_corner_S_neumann_E_dirichlet;
    grid_p.is_corner_N_neumann_E_dirichlet;
    grid_p.is_concave_NWS;                      % 凹角分类
    grid_p.is_concave_NES;
    grid_p.is_concave_WNE;
    grid_p.is_concave_WSE;
};

% 测试1：覆盖检查 - 每个节点都必须被分类
total_mask = false(Ny, Nx);
for i = 1:length(all_masks_for_check)
    total_mask = total_mask | all_masks_for_check{i};
end
if all(total_mask(:))
    coverage_ok = true;
    fprintf('[通过] 覆盖检查: 所有网格节点都已分类。\n');
else
    unclassified_count = nnz(~total_mask);
    coverage_ok = false;
    fprintf('[失败] 覆盖检查: 发现 %d 个未分类节点！\n', unclassified_count);
end

% 测试2：唯一性检查 - 没有重叠分类
overlap_map = zeros(Ny, Nx, 'uint8');
for i = 1:length(all_masks_for_check)
    overlap_map = overlap_map + uint8(all_masks_for_check{i});
end
if all(overlap_map(:) <= 1)
    uniqueness_ok = true;
    fprintf('[通过] 唯一性检查: 所有节点分类都是唯一的。\n\n');
else
    overlapped_count = nnz(overlap_map > 1);
    uniqueness_ok = false;
    fprintf('[失败] 唯一性检查: 发现 %d 个节点有重叠分类！\n\n', overlapped_count);
end

%% 用户确认界面
set(fig_grid, 'MenuBar', 'figure', 'ToolBar', 'figure');
uicontrol(fig_grid, 'Style', 'text', 'String', sprintf('网格: %d×%d | 网格尺寸: %.4f m | 计算域: %.3f×%.3f m', Nx, Ny, h, L, H), 'Position', [50, 50, 600, 25], 'FontSize', 11, 'HorizontalAlignment', 'left', 'BackgroundColor', get(fig_grid, 'Color'));
if coverage_ok && uniqueness_ok, status_str = '✅ 完整性检查: 通过'; status_color = [0.2, 0.6, 0.2];
else, status_str = '⚠️ 完整性检查: 失败 (详见命令窗口)'; status_color = [0.8, 0.1, 0.1]; end
uicontrol(fig_grid, 'Style', 'text', 'String', status_str, 'Position', [50, 20, 400, 25], 'FontSize', 11, 'FontWeight', 'bold', 'ForegroundColor', status_color, 'BackgroundColor', get(fig_grid, 'Color'));
uicontrol(fig_grid, 'Style', 'pushbutton', 'String', '✅ 继续仿真', 'Position', [fig_grid.Position(3)-370, 25, 180, 45], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.2, 0.7, 0.2], 'ForegroundColor', 'white', 'Callback', @continue_callback);
uicontrol(fig_grid, 'Style', 'pushbutton', 'String', '🔄 返回设置', 'Position', [fig_grid.Position(3)-170, 25, 150, 45], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.7, 0.3, 0.1], 'ForegroundColor', 'white', 'Callback', @retry_callback);

    % 用户交互的回调函数
    function continue_callback(~, ~), continue_flag = true; close(fig_grid); end
    function retry_callback(~, ~), continue_flag = false; close(fig_grid); end
    function close_callback(~, ~)
        if isempty(continue_flag), continue_flag = false; end
        delete(fig_grid);
    end
set(fig_grid, 'CloseRequestFcn', @close_callback);
uiwait(fig_grid);                               % 等待用户决定

end
