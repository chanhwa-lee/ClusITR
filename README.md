# Individualized Treatment Rule under Clustered Interference

## Replication code for Zhang & Imai (2023)[https://arxiv.org/abs/2311.02467]

- Mixed Integer Linear Programming to find optimal ITR under Clustered Interference
- Two plots are
  1. Policy Value $$V(\hat{\pi})$$ for each learned policy $$\hat{\pi}$$ by different methods -- see how the regret is decreasing
  2. Learned $$\hat{\beta}$$ values compared to the Oracle $\beta$ values -- see how the learned policy coefficients are close to the oracle coefficients

## Future direction

- Extension to survival setting (right censored outcome) + observational study (unknown propensity)
- Find DGP with semiparametric heterogeneous additive outcome model assumption
  
$$
E[T_{ij}|A_i,X_i] = g_j^{(0)}(X_i) + \sum_k g_j^{(k)}(X_i) A_{ik}
$$

- For this, $g_j^{(k)}(X_i)$ should be positive










## Reference

Zhang, Yi, and Kosuke Imai. "Individualized Policy Evaluation and Learning under Clustered Network Interference." arXiv preprint arXiv:2311.02467 (2023).
