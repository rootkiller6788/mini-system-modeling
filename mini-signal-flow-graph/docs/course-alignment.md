# Course Alignment — mini-signal-flow-graph

## Nine-School Curriculum Mapping

| School | Course | SFG Coverage | Our Module |
|--------|--------|-------------|------------|
| **MIT** | 6.302 Feedback Systems | Mason's rule, block diagram ↔ SFG | `sfg_mason.h/c` |
| **MIT** | 6.241 Dynamic Systems | State-space ↔ SFG, canonical forms | `sfg_state_space.h/c` |
| **Stanford** | ENGR105 Feedback Control | SFG reduction, TF computation | `sfg_reduction.h/c` |
| **Stanford** | ENGR207B Linear Systems | MIMO SFG, transfer function matrix | `sfg_mason.h/c` (MIMO) |
| **Berkeley** | ME132 Dynamic Systems | Mechanical system SFG modeling | `example_mass_spring.c` |
| **Caltech** | CDS 101/110 Intro | SFG basics, Mason's rule introduction | `sfg_core.h/c` |
| **ETH** | 151-0591 Control I | SFG reduction rules R1-R5 | `sfg_reduction.h/c` |
| **Cambridge** | 3F2 Systems & Control | Sensitivity, stability via SFG | `sfg_analysis.h/c` |
| **Georgia Tech** | AE 6530 Optimal | MIMO analysis, sensitivity matrix | `sfg_analysis.h/c` |
| **Purdue** | ECE 602 Lumped Systems | Circuit analysis via SFG | `example_electrical.c` |
| **清华** | 自动控制原理 | Mason公式, SFG化简, 框图转换 | All modules |

## Textbook Alignment

| Textbook | Chapter | Our Implementation |
|----------|---------|-------------------|
| Ogata, Modern Control Engineering (5th Ed.) | Ch.3 (SFG), Ch.9 (State-Space) | `sfg_mason`, `sfg_state_space` |
| D'Azzo & Houpis, Linear Control System (4th Ed.) | Ch.3 (SFG & Mason) | `sfg_mason`, `sfg_reduction` |
| Kuo, Automatic Control Systems (9th Ed.) | Ch.5 (State-Variable) | `sfg_state_space` |
| Mason (1953, 1956) IRE Papers | Original SFG theory | `sfg_mason` (algorithm) |
| Johnson (1975) SIAM | Cycle enumeration | `sfg_path` (loop finding) |
| Lorens (1964) Flowgraphs | SFG reduction rules | `sfg_reduction` |
