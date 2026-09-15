# behavioral-hetero
This is the repository for the NSF IHBEM behavioral heterogeneity project. This repository contains the data presented in and the code used for all simulations and analysis in the publication [1].

All code was run using MATLAB R2025b.

Additional instructions:
- Figures from the data will generate into the empty "outputs_data" folder upon running "run_load_and_plot_data.m"
- Figures for the base two-group model (Figures 3-7, A1, B1, C1-C6, and D1) will generate into the empty "outputs_seirb_twogroup" folder upon running "run_seirb_twogroup.m"
- Figures for the multi-group model (Figures 8 and 9) will generate into the empty "outputs_seirb_multigroup" folder upon running "run_seirb_multigroup.m"
- Data on relative risk by age group is contained in "COVID-19_risk_by_age.xlsx" and used by file "run_load_and_plot_data.m"
- Data on the contact matrix is contained in the "Excel_1_(USA).xlsx" file and used by file "run_seirb_multigroup.m"

[1] Fattahpour H, Ghaffarzadegan N, Childs LM. Not everyone panics alike: When heterogeneous behavioral responses change epidemic outcomes. Journal of Theoretical Biology. doi: 10.1016/j.jtbi.2026.112590
