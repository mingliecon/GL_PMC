# Replication Package for: Identification of Semiparametric Panel Multinomial Choice Models with Infinite-Dimensional Fixed Effects

Wayne Yuan Gao and Ming Li  
June 10, 2025

This replication package accompanies Gao and Li, *“Identification of Semiparametric Panel Multinomial Choice Models with Infinite-Dimensional Fixed Effects.”*

Empirical datasets are not uploaded on this repository due to confidentiality
---

## 1. Description of package content

The replication package includes all raw data and codes for replicating the tables and figures in the paper.  

- Two main folders:  
  - `/simulation` – contains files for replicating Tables 1–2, Tables 7–8, Figures 2–3.  
  - `/empirical` – contains files for replicating Tables 3–6.  

Each folder has a master code:  
- `run_simulation.m` in `/simulation` generates all simulation tables and figures (saved in `/simulation/result_simulation`).  
- `run_empirical.m` in `/empirical` generates all empirical tables (saved in `/empirical/result_empirical`).  

Subfolders are also provided if you want to generate specific results:  
- `/simulation/table_1_2`, `/simulation/table_7`, `/simulation/table_8`, `/simulation/figure_2_3`  
- `/empirical/table_3`, `/empirical/table_4_5`, `/empirical/table_6`  

Inside each subfolder:  
- `src_*` folder contains the execution code (e.g., `run_table_1_2.m`).  
- `result_*` folder contains the outputs.  

Example: Running `run_table_1_2.m` in `/simulation/table_1_2/src_table_1_2` generates Tables 1 and 2 inside `/simulation/table_1_2/result_table_1_2` and `/simulation/result_simulation`.

---

## 2. Code execution instructions

**Simulation results**  
1. Run `run_simulation.m` in `/simulation`.  
2. Outputs: Tables 1–2, Tables 7–8, Figures 2–3 → `/simulation/result_simulation`.  

**Empirical results**  
1. Run `run_empirical.m` in `/empirical`.  
2. Outputs: Tables 3–6 → `/empirical/result_empirical`.  

**Specific tables or figures**  
- Run the corresponding `run_*.m` file in the respective subfolder.  
- Results will be saved in the associated `result_*` folder.  

**Parallel Processing**  
- Run with multiple cores using: `run_empirical(N)` (where `N` = number of cores).  
- Example: `run_empirical(12)` for 12 cores.  
- Check available cores in MATLAB with `feature('numcores')` or `maxNumCompThreads`.

---

## 3. Computational requirements

**Hardware used:**  
- (Specify your hardware setup if needed)

**Software:**  
- MATLAB Version 24.2 (R2024b)  
- Curve Fitting Toolbox Version 24.2  
- Optimization Toolbox Version 24.2  
- Parallel Computing Toolbox Version 24.2  
- Statistics and Machine Learning Toolbox Version 24.2  

---

## 4. List of tables and figures

### Simulation Results
| Table/Figure | Code | Folder | Output Folder | Notes |
|--------------|------|--------|---------------|-------|
| Table 1–2 | `run_simulation.m` | `/simulation/table_1_2` | `/simulation/table_1_2/result_table_1_2` and `/simulation/result_simulation` | Results split into upper and bottom parts |
| Table 7 | `run_simulation.m` | `/simulation/table_7` | `/simulation/table_7/result_table_7` and `/simulation/result_simulation` |  |
| Table 8 | `run_simulation.m` | `/simulation/table_8` | `/simulation/table_8/result_table_8` and `/simulation/result_simulation` |  |
| Figures 2–3 | `run_simulation.m` | `/simulation/figure_2_3` | `/simulation/figure_2_3/result_figure_2_3` and `/simulation/result_simulation` |  |

### Empirical Results
| Table/Figure | Code | Folder | Output Folder |
|--------------|------|--------|---------------|
| Table 3 | `run_empirical.m` | `/empirical/table_3` | `/empirical/table_3/result_table_3` and `/empirical/result_empirical` |
| Table 4–5 | `run_empirical.m` | `/empirical/table_4_5` | `/empirical/table_4_5/result_table_4_5` and `/empirical/result_empirical` |
| Table 6 | `run_empirical.m` | `/empirical/table_6` | `/empirical/table_6/result_table_6` and `/empirical/result_empirical` |

---
