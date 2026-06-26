# Course Alignment — mini-fluid-thermal-modeling

## Nine-School Curriculum Mapping

### MIT
| Course | Topic | This Module |
|--------|-------|-------------|
| 2.005 Thermal-Fluids Engineering I | Bernoulli, pipe flow, heat exchangers | `fluid_mechanics.c`, `heat_transfer.c` |
| 2.006 Thermal-Fluids Engineering II | Convection, radiation, HX design | `heat_transfer.c` Nusselt correlations |
| 2.25 Advanced Fluid Mechanics | Navier-Stokes, boundary layers | `fluid_mechanics.c` analytical solutions |
| 2.51 Intermediate Heat Transfer | Fin analysis, transient conduction | `heat_transfer.c` sections |

### Stanford
| Course | Topic | This Module |
|--------|-------|-------------|
| ME 346A Heat Transfer | Conduction, convection, radiation fundamentals | `heat_transfer.c` |
| ME 351A Fluid Mechanics | Pipe flow, dimensional analysis | `fluid_mechanics.c` |
| ME 352C Advanced Heat Transfer | ε-NTU, compact HX | `heat_transfer.c` HX methods |

### UC Berkeley
| Course | Topic | This Module |
|--------|-------|-------------|
| ME 105 Thermodynamics | Property relations, cycles | `fluid_thermal_core.c` fluid properties |
| ME 106 Fluid Mechanics | Bernoulli, NS, Moody | `fluid_mechanics.c` |
| ME 109 Heat Transfer | Conduction, fins, HX | `heat_transfer.c` |

### Michigan
| Course | Topic | This Module |
|--------|-------|-------------|
| ME 320 Fluid Mechanics | Laminar/turbulent flow, pipe systems | `fluid_mechanics.c` |
| ME 420 Heat Transfer | Natural/forced convection | `heat_transfer.c` correlations |
| AERO 533 Combustion | Compressible flow fundamentals | `fluid_thermal_core.c` Mach/compressibility |

### Purdue
| Course | Topic | This Module |
|--------|-------|-------------|
| ME 505 Intermediate Heat Transfer | Multi-mode, fins, HX methods | `heat_transfer.c` |
| ME 509 Intermediate Fluid Mechanics | BL theory, potential flow | `fluid_mechanics.c` BL section |
| ME 597 Plasma Engineering | Stefan-Boltzmann, radiation | `heat_transfer.c` radiation section |

### TU Munich
| Course | Topic | This Module |
|--------|-------|-------------|
| MW 0798 Thermodynamics | Fluid properties, cycles | `fluid_thermal_core.c` |
| MW 0854 Heat and Mass Transfer | Conduction, convection, radiation | `heat_transfer.c` |
| MW 0860 Fluid Mechanics | Pipe flow, pumps, hydraulic systems | `thermal_fluid_systems.c` |

### ETH Zurich
| Course | Topic | This Module |
|--------|-------|-------------|
| 151-0103 Fluid Dynamics | NS equations, exact solutions | `fluid_mechanics.c` NS section |
| 151-0111 Heat Transfer | Full heat transfer spectrum | `heat_transfer.c` |
| 227-0455 Engineering EM | Stefan-Boltzmann, EM waves | `heat_transfer.c` radiation |

### Tsinghua University
| Course | Topic | This Module |
|--------|-------|-------------|
| 工程热力学 | Property tables, cycles | `fluid_thermal_core.c` |
| 传热学 | Full conduction/convection/radiation | `heat_transfer.c` |
| 流体力学 | Pipe flow, BL, dimensional analysis | `fluid_mechanics.c` |

### Cambridge
| Course | Topic | This Module |
|--------|-------|-------------|
| 4A3 Fluid Mechanics | NS, BL theory, turbulence intro | `fluid_mechanics.c` |
| 4A9 Heat Transfer | Conduction, convection, exchangers | `heat_transfer.c` |
| 4M12 Thermal Systems | Combined thermal-fluid systems | `thermal_fluid_systems.c` |

### Oxford
| Course | Topic | This Module |
|--------|-------|-------------|
| B1 Fluid Mechanics | Laminar flow, pipe networks | `fluid_mechanics.c` + `numerical_fluid_thermal.c` |
| B2 Heat Transfer | Conduction, fins, transient | `heat_transfer.c` |

### Caltech
| Course | Topic | This Module |
|--------|-------|-------------|
| Ae/APh/CE/ME 101a Fluid Mechanics | NS exact solutions, BL | `fluid_mechanics.c` |
| Ae/ME 103 Heat Transfer | Combined modes, HX | `heat_transfer.c` |
