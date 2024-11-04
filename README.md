# Individualized Treatment Rule under Clustered Interference

Overleaf document: https://www.overleaf.com/5147233984szhkrjpwdfrb#6bcbe7

## Replication code for Zhang & Imai (2023)[https://arxiv.org/abs/2311.02467]

- See `simul_Zhang_replicate.R`: simulation results are in `Output_Zhang`

- Mixed Integer Linear Programming to find optimal ITR under Clustered Interference
![MIP formulation Note](https://github.com/chanhwa-lee/ClusITR/blob/main/MIP%20formulation%20Note.jpeg)

- Two plots will be generated if running last part of `simul_Zhang_replicate.R`:
  1. Policy Value $$V(\hat{\pi})$$ for each learned policy $$\hat{\pi}$$ by different methods -- see how the regret is decreasing

   ![Learned Policy Value by Method](https://github.com/user-attachments/assets/0c6dfec4-6052-4a0d-af8b-172c79c61602)

     
  2. Learned $$\hat{\beta}$$ values compared to the Oracle $\beta$ values -- see how the learned policy coefficients are close to the oracle coefficients
     
 ![Bias of Beta Estimates by Method and Sample Size](https://github.com/user-attachments/assets/57ac7af4-3baf-4d0d-bd8d-3e6f357950e2)

      

## Future direction

- Extension to survival setting (right censored outcome) + observational study (unknown propensity)
- Find DGP with semiparametric heterogeneous additive outcome model assumption
  
$$
E[T_{ij}|A_i,X_i] = g_j^{(0)}(X_i) + \sum_k g_j^{(k)}(X_i) A_{ik}
$$

- For this, $g_j^{(k)}(X_i)$ should be positive










## Reference

Zhang, Yi, and Kosuke Imai. "Individualized Policy Evaluation and Learning under Clustered Network Interference." arXiv preprint arXiv:2311.02467 (2023).
