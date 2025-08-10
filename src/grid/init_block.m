function [block_info, L, h, Nx, Ny] = init_block(image_filename, H, node_scale)
%INIT_BLOCK 从图像中提取障碍物几何形状并生成网格参数
%   处理二值障碍物图像以提取边界坐标并建立计算网格。
%   支持多个障碍物，具有自动边界跟踪和坐标变换功能。
%
%   输入参数:
%       image_filename - 障碍物图像文件路径(PNG, JPG等)
%       H              - 物理计算域高度
%       node_scale     - 目标网格节点总数(~Nx*Ny)
%
%   输出参数:
%       block_info - 障碍物结构的元胞数组，每个包含:
%                    .points   - [N×2]数组的边界坐标
%                    .image    - 用于可视化的原始图像切片
%                    .mask     - 用于几何查询的二值掩码
%                    .x_coords - 物理x范围[xmin, xmax] 
%                    .y_coords - 物理y范围[ymin, ymax] 
%       L          - 物理计算域宽度[m](根据图像纵横比)
%       h          - 网格间距[m](两个方向均匀)
%       Nx, Ny     - x和y方向的网格单元数
%
%   图像处理:
%   - 将彩色图像转换为二值图像(阈值=128)
%   - 应用形态学闭运算填补像素间隙
%   - 移除小于10像素的噪声颗粒
%   - 仅跟踪外边界(无内部孔洞)
%
%   See also BWBOUNDARIES, IMCLOSE, CONFIRM_GRID.

%% 图像加载和预处理
try
    if ~exist(image_filename, 'file')
        error('文件未找到: %s', image_filename);
    end 
    imageMatrix = imread(image_filename);
catch ME
    error('加载图像"%s"失败.\n原因: %s\n请检查文件路径和格式(png, jpg等).', ...
          image_filename, ME.message);
end

% 转换为灰度图像以便一致处理
if size(imageMatrix, 3) == 3
    img_gray = rgb2gray(imageMatrix);            % RGB到灰度转换
else
    img_gray = imageMatrix;                      % 已经是灰度图
end

% 二值化和噪声清理
img_bw = img_gray < 128;                         % 阈值二值化(黑色=障碍物)

% 形态学操作改善边界质量
se = strel('disk', 1);                           
img_bw = imclose(img_bw, se);                    % 闭合1像素间隙
img_bw = bwareaopen(img_bw, 10);                 % 移除小噪声颗粒

% 根据图像纵横比计算物理计算域宽度
[imgHeight, imgWidth] = size(img_bw);
L = H * (imgWidth / imgHeight);                  % 保持图像纵横比

%% 障碍物边界提取
% 仅提取外边界(忽略内部孔洞)
[B, ~] = bwboundaries(img_bw, 'noholes');

% 处理未检测到障碍物的情况
if isempty(B)
    block_info = {};
    fprintf('警告: 在图像"%s"中未找到障碍物.\n', image_filename);
    return;
end

%% 处理每个检测到的障碍物
numObstacles = length(B);
block_info = cell(numObstacles, 1);

for k = 1:numObstacles
    boundary_pixels = B{k};                      % 像素坐标中的边界
    
    % 将像素坐标转换为物理坐标
    pixel_x = boundary_pixels(:, 2);             % 列索引(x方向)
    pixel_y = boundary_pixels(:, 1);             % 行索引(y方向)
    scaled_x = (pixel_x / imgWidth) * L;         % 缩放到物理x坐标
    scaled_y = H - (pixel_y / imgHeight) * H;    % 缩放并翻转y轴(图像约定)
    scaled_points = [scaled_x, scaled_y];
    
    % 为可视化和几何查询提取边界框
    min_row = min(pixel_y); max_row = max(pixel_y);
    min_col = min(pixel_x); max_col = max(pixel_x);
    block_image_slice = imageMatrix(min_row:max_row, min_col:max_col, :);
    block_mask_slice = img_bw(min_row:max_row, min_col:max_col);
       
    % 计算边界框的物理坐标范围
    x_coords = [min_col/imgWidth, max_col/imgWidth] * L;
    y_coords = H - [max_row/imgHeight, min_row/imgHeight] * H;
    
    % 以结构化格式存储障碍物数据
    block_info{k} = struct('points', scaled_points, ...
                           'image', block_image_slice, ...
                           'mask', block_mask_slice, ...
                           'x_coords', x_coords, ...
                           'y_coords', y_coords);
end

%% 生成均匀网格参数
L = H * (imgWidth / imgHeight);                  % 计算域宽度
h = sqrt((L * H) / node_scale);                  % 均匀网格间距
Nx = round(L / h);                               % x方向网格单元
Ny = round(H / h);                               % y方向网格单元

end
