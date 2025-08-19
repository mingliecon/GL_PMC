# Repository of Replication Codes for: Identification of Semiparametric Panel Multinomial Choice Models with Infinite-Dimensional Fixed Effects

**Authors:** Wayne Yuan Gao and Ming Li  
**Last Updated:** August 2025  

This repository provides the replication codes for Gao and Li, *“Identification of Semiparametric Panel Multinomial Choice Models with Infinite-Dimensional Fixed Effects.”*  

It contains all scripts and output files necessary to reproduce the tables and figures reported in the paper.  
> ⚠️ **Note**: Due to confidentiality restrictions, the empirical datasets are not included in this repository. As a result, Tables 3–5 cannot be fully replicated, though the corresponding code is provided.


---

## Content Page of README

1. [Introduction of the Project](#1-introduction-of-the-project)  
2. [Description of Package Content](#2-description-of-package-content)  
3. [Computational Requirements](#3-computational-requirements)  

---

## 1. Introduction of the project
we could add a description of the project, such as the abstract of the paper here ???

---


## 2. Description of package content

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

```plaintext
Replication-Package/
│
├── algorithm/
│
├── paper/
│
├── empirical/
│   ├── table_4_5/
│   │   ├── src_table_4_5/
│   │   │   ├──create_interaction_variables.m
│   │   │   ├──Qfunc.m
│   │   │   ├──Table_4_nlogn.m
│   │   │   ├──Table_4.m
│   │   │   └──run_table_4_5.m
│   │   │ 
│   │   └── result_table_4_5/
│   │       ├──table4.csv
│   │       └──table5.csv
│   │
│   ├── table_6/
│   │   ├── src_table_6/
│   │   │   ├──GL
│   │   │   │   ├──create_interaction_variables.m
│   │   │   │   ├──Qfunc.m
│   │   │   │   ├──Gl_Table6.m
│   │   │   │   └──run_GL.m
│   │   │   └──run_table_6.m
│   │   │  
│   │   └── result_table_6/
│   │       └──table6.csv
│   │
│   ├── result_empirical/
│   │   ├──table4.csv
│   │   ├──table5.csv
│   │   └──table6.csv
│   │
│   └── run_empirical.m
│
├── simulation/
│   ├── table_1_2/
│   │   ├── src_table_1_2/
│   │   │   ├──create_interaction_variables.m
│   │   │   ├──main_mc_agg.m
│   │   │   ├──mc_agg_0_setup.m
│   │   │   ├──mc_agg_1_dgp.m
│   │   │   ├──mc_agg_2_gamma.m
│   │   │   ├──mc_agg_3_beta_D3.m
│   │   │   ├──mc_agg_4_evaluation.m
│   │   │   ├──Qfunc.m
│   │   │   └──run_table_1_2.m
│   │   │
│   │   └── result_table_1_2/
│   │       ├──table1.csv
│   │       ├──table2_bottom_part.csv
│   │       └──table2_upper_part.csv
│   │
│   ├── table_7/
│   │   ├── src_table_7/
│   │   │   ├──PointNID
│   │   │   │     ├──create_interaction_variables.m
│   │   │   │     ├──run_pointnid.m
│   │   │   │     ├──Qfunc.m
│   │   │   │     └──simulation_nid.m
│   │   │   │ 
│   │   │   ├──PointID
│   │   │   │     ├──create_interaction_variables.m
│   │   │   │     ├──run_pointid.m
│   │   │   │     ├──Qfunc.m
│   │   │   │     └──simulation_id.m
│   │   │   │
│   │   │   └──run_table_7.m
│   │   │
│   │   └── result_table_7/
│   │       └──table7.csv
│   │
│   ├── table_8/
│   │   ├── src_table_8/
│   │   │   ├──create_interaction_variables.m
│   │   │   ├──MC_AGG_20230909.m
│   │   │   ├──Qfunc.m
│   │   │   └──run_table_8.m
│   │   │
│   │   └── result_table_8/
│   │       └──table8.csv
│   │
│   ├── figure_2_3/
│   │   ├── src_figure_2_3/
│   │   │   ├──code_fig_2_3.m
│   │   │   ├──run_fig_2_3.m
│   │   │   └──Qfunc.m
│   │   │ 
│   │   └── result_figure_2_3/
│   │       ├──Sim_3D_Grid_sieve_theta_1.png
│   │       └──Sim_3D_trues_beta_1.png
│   │
│   ├── result_simulation/
│   │   ├──table1.csv
│   │   ├──table2_bottom_part.csv
│   │   ├──table2_upper_part.csv
│   │   ├──table7.csv
│   │   └──table8.csv
│   │   
│   └── run_simulation.m
│
└── README.md
```

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

