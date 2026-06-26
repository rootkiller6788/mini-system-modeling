# Mini System Modeling（迷你系统建模）

**从零开始、零依赖的 C 语言实现**，涵盖大学级系统建模与控制理论。每个模块对应 MIT（及其他顶尖大学）的一门或多门课程，将教科书中的公式转化为可运行的 C 代码，实现理论与实践的桥接。

## 子模块

| 子模块 | 主题 | 参考课程 |
|--------|------|----------|
| [mini-block-diagram](mini-block-diagram/) | 框图数据结构、框图代数、系统化化简、梅森增益公式、信号流图 | MIT 6.302, Stanford ENGR105, ETH 151-0591 |
| [mini-electrical-modeling](mini-electrical-modeling/) | 电路拓扑（KVL/KCL）、RLC 元件、状态空间模型、机电系统、双端口网络 | MIT 6.302, Stanford ENGR105, Berkeley EE 219 |
| [mini-fluid-thermal-modeling](mini-fluid-thermal-modeling/) | 流体力学（伯努利、纳维-斯托克斯）、传热（传导、对流、辐射）、数值方法、耦合热流系统 | MIT 2.006, MIT 2.51, Stanford ME 351B |
| [mini-linearization](mini-linearization/) | 泰勒展开、雅可比线性化、描述函数法、反馈线性化、增益调度、Koopman 算子 | MIT 6.241J, MIT 2.154, Stanford ENGR207B |
| [mini-mechanical-modeling](mini-mechanical-modeling/) | 质量-弹簧-阻尼元件、传递函数、状态空间、多自由度模态分析、旋转系统、悬架与隔振 | MIT 2.003, MIT 2.004, Stanford ENGR105 |
| [mini-mechatronic-modeling](mini-mechatronic-modeling/) | 执行器选型、机电耦合、传感器建模、机电状态空间、集成系统设计 | MIT 2.154, Stanford ENGR105, ETH 151-0591 |
| [mini-signal-flow-graph](mini-signal-flow-graph/) | 信号流图核心结构、梅森增益公式、路径/回路枚举、拓扑化简、矩阵方法、状态空间转换、灵敏度分析 | MIT 6.302, Stanford ENGR105, Cambridge 3F2 |
| [mini-transfer-function](mini-transfer-function/) | 多项式代数、s 域有理函数、互联代数、Routh-Hurwitz 判据、Nyquist、Bode、SS↔TF 转换、模型降阶 | MIT 6.302, Stanford ENGR105, Berkeley ME132, ETH 151-0591 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有课程对齐说明
- **实用演示程序** — 框图化简、电路仿真、流体热分析、线性化工作流、机械系统辨识、机电设计、信号流图分析、传递函数评估

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-block-diagram
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-system-modeling/
├── mini-block-diagram/           # 框图代数、化简、梅森公式
├── mini-electrical-modeling/     # 电路拓扑、RLC 元件、双端口网络
├── mini-fluid-thermal-modeling/  # 流体力学、传热、耦合系统
├── mini-linearization/           # 泰勒级数、雅可比、描述函数、增益调度
├── mini-mechanical-modeling/     # 质量-弹簧-阻尼、多自由度、旋转系统
├── mini-mechatronic-modeling/    # 执行器、传感器、机电耦合、系统集成
├── mini-signal-flow-graph/       # 信号流图核心、梅森公式、路径枚举、化简
└── mini-transfer-function/       # 多项式代数、稳定性分析、模型转换、降阶
```

## 许可证

MIT
