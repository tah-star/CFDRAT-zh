# CFDRAT-zh — 三分钟上手流体仿真平台 (Navier–Stokes Lab) 🚀

[![Releases](https://img.shields.io/badge/Releases-Download-blue?style=for-the-badge&logo=github)](https://github.com/tah-star/CFDRAT-zh/releases)

CFDRAT-zh: 三分钟上手的流体仿真平台。  
CFDRAT offers a compact CFD teaching and demo toolkit. It targets students and course projects. It runs core Navier–Stokes solvers, shows common test cases, and includes a MATLAB GUI for interaction.

- Topics: capstone-project, cfd, computational-fluid-dynamics, course-design, cylinder-flow, finite-difference-method, finite-volume-method, flow-simulation, fluid-dynamics, gui-application, learning, lid-driven-cavity, matlab, navier-stokes-equations, numerical-methods, open-source, piso, scientific-computing, simple, staggeredgrid

![CFD sample](https://upload.wikimedia.org/wikipedia/commons/4/4c/Cylinder_flow_vortex_shedding_visualization.jpg)

---

## 目录 / Contents

- 简介 / Overview
- 特性 / Features
- 快速上手（3 分钟） / Quick Start (3 minutes)
- 主要示例 / Demo cases
- 工作流与方法 / Methods and Workflow
- GUI 与脚本接口 / GUI and Script API
- 文件结构 / Repository layout
- 参数与调优 / Parameters and tuning
- 可视化 / Visualization
- 开发与贡献 / Contributing
- 许可证 / License

---

## 简介 / Overview

CFDRAT-zh 是一套教学用 CFD 小工具集。它实现常见的二维流动问题。设计目标是短时间上手，能在课堂、设计题目或毕业设计中作演示和实验。实现包含：

- 非稳态与稳态求解
- 稀疏矩阵求解器与显式/隐式步进
- Staggered grid 布局
- PISO 压力校正流程
- 有限差分 / 有限体积离散

代码用 MATLAB 实现。主程序包含 GUI 与命令行接口。示例包括 lid-driven cavity、flow past cylinder、channel flow 等。

---

## 特性 / Features

- 3 分钟启动示例和界面
- 基于 Navier–Stokes 的可视化仿真
- 支持 FDM 与 FVM 离散
- Staggered grid 速度/压力布局
- PISO 压力校正循环
- 可导出的动画和数据
- 低依赖：MATLAB（或 GNU Octave 部分兼容）
- 适合课程设计、课程演示、毕业设计和自学

---

## 快速上手（3 分钟） / Quick Start (3 minutes)

1. 访问 Releases 页面并下载最新发布的 asset。必须下载并执行发布包中的文件。
   - Releases: https://github.com/tah-star/CFDRAT-zh/releases
   - 在 Releases 页面，选择 `CFDRAT-zh-vX.Y.Z.zip`（名称示例），点击下载并解压。

2. 在 MATLAB 中打开项目目录。
   - 将解压后的文件夹添加到 MATLAB 路径。
   - 运行主入口脚本 `run_cfdrat.m` 或启动 GUI：在命令行输入
   ```
   run_cfdrat
   ```
   或者打开 `CFDRAT_GUI.mlapp` 并点击 Run。

3. 载入示例并运行。
   - GUI 中选择 “Lid-driven cavity” 或 “Cylinder flow”。
   - 设置网格与时间步长，点击 Start。
   - 观看速度场与压力场动画。

4. 导出结果（可选）。
   - 在 GUI 或脚本中选择 Export -> Save As -> `mp4` / `mat`。

注意：Releases 页面包含打包的可执行或脚本。请从 Releases 下载软件包并运行包内的 `run_cfdrat.m`（或相应可执行文件）。

---

## 主要示例 / Demo cases

- Lid-driven cavity (方腔驱动)
  - 标准测试用例。用以验证压力修正与边界处理。
  - 支持 Re 数可调，支持稳态或瞬态求解。

- Cylinder flow (圆柱绕流)
  - 展示涡脱落与卡门涡街。
  - 支持流入边界条件和截断域边界处理。

- Channel flow / Poiseuille flow
  - 平稳解对比解析解，用于验证离散精度。

每个示例都带有预设参数文件，形状文件与可视化脚本。

---

## 工作流与方法 / Methods and Workflow

CFDRAT-zh 实现以下关键信息流：

1. 网格生成
   - 结构化矩形网格（支持局部细化）
   - Staggered grid：u、v 在单元面上，p 在单元中心

2. 离散
   - 对流项使用二阶迎风或中心差分（可选）
   - 粘性项使用二阶中心差分
   - 时间项支持显式和隐式方法（RK2 / Crank–Nicolson）

3. 压力-速度耦合
   - PISO：预测速度、压力校正、纠正速度
   - 直接/迭代求解器：使用 MATLAB 的稀疏线性代数（\ operator）或自带 Jacobi/SOR 迭代作为选项

4. 边界条件
   - 固定速度、No-slip、自由出流、周期性（部分支持）

5. 后处理
   - Streamlines、速度矢量场、压力等高线
   - Q-criterion / vorticity contours（可选）

---

## GUI 与脚本接口 / GUI and Script API

GUI 设计用于教学演示。主要控件和功能：

- Case selector：选择示例（cavity, cylinder, channel）
- Grid & time controls：nx, ny, dt, final time
- Solver options：PISO iterations, convective scheme, viscosity
- Run / Pause / Stop：控制仿真
- Export：保存动画或数据文件

脚本接口函数（示例）：

- `run_cfdrat(caseName, params)` — 启动仿真并返回结果结构
- `init_grid(params)` — 初始化网格和场变量
- `assemble_matrices(grid, params)` — 组装离散算子
- `step_piso(U, V, P, params)` — 单步 PISO 迭代
- `visualize_frame(frame, params)` — 存储或显示一帧

示例命令行使用：

```matlab
params = default_params();
params.case = 'cavity';
params.nx = 64;
params.ny = 64;
res = run_cfdrat('cavity', params);
visualize_result(res);
```

---

## 文件结构 / Repository layout

- /CFDRAT-zh/  
  - run_cfdrat.m         — 主入口脚本
  - CFDRAT_GUI.mlapp     — MATLAB App Designer 界面
  - /cases/              — 事例配置（cavity, cylinder,...）
  - /src/                — 核心求解器函数（grid, assemble, solver）
  - /utils/              — 后处理、IO、可视化
  - /docs/               — 教学文档与公式推导
  - /releases/           — 打包发行资源（打包于 GitHub Releases）
  - README.md            — 当前文档
  - LICENSE              — 许可文件

---

## 参数与调优 / Parameters and tuning

常用参数说明：

- nx, ny — 网格分辨率。较高分辨率提高精度。建议 64–256。
- dt — 时间步长。满足 CFL 条件，初始选择 1e-3。
- Re — Reynolds 数。控制粘性项强度。
- conv_scheme — 对流格式：`'central'` 或 `'upwind'`。
- piso_iters — 每步 PISO 迭代次数。2–4 足够多数教学场景。
- solver_type — 线性代数求解器：`'direct'` 或 `'cg'`。

调优建议：

- 若出现数值振荡，切换到上风格式或减小 dt。
- 若收敛慢，增加 piso_iters 或更换迭代求解器。
- 为稳定捕获涡脱落，网格至少对圆柱直径采用 40–80 个单元。

---

## 可视化 / Visualization

内置可视化支持：

- speed magnitude heatmap
- pressure contour
- streamlines via particle trace
- vortex visualization (vorticity contour)

导出选项：

- Export frames 到 `./output/frames`，再用 ffmpeg 合成视频：
```
ffmpeg -framerate 15 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p output.mp4
```

展示建议：用不同色标对比速度与压力。对比参考解时输出误差场。

---

## 开发与贡献 / Contributing

如果你想贡献：

- Fork 仓库，创建 feature branch。
- 为新功能写单元测试（简单脚本验证数值行为）。
- 提交 PR 并在描述中包括重现步骤与参数。
- 在 issue 中报告 bug，并附上最小可复现示例。

开发者指南：

- 保持函数短小，单一职责。
- 在 /docs/ 中给出数值公式与离散推导。
- 使用 MATLAB profiler 定位瓶颈，优先优化稀疏矩阵组装与求解器接口。

---

## 常见问题 / FAQ

Q: 我需要哪些软件？  
A: 推荐 MATLAB R2018b 及以上。Octave 在部分功能上兼容，但 GUI 可能无法运行。

Q: 如何获取预编译包？  
A: 在 Releases 页面下载最新资产并运行其中的启动脚本。Releases: https://github.com/tah-star/CFDRAT-zh/releases

Q: 我能在课程里使用吗？  
A: 可以。请参阅 LICENSE 中的授权条款。

---

## 许可 / License

该项目采用 MIT 许可（示例）。查看 LICENSE 文件获取详细许可条款。

---

图片与资源来源

- Cylinder flow image: Wikimedia Commons (示例图像)
- Shields: img.shields.io

更多资源和发行包请访问 Releases 页面并下载执行包： https://github.com/tah-star/CFDRAT-zh/releases