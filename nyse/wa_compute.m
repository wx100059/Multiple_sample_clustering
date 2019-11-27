function wasserstein_distance = wa_compute(mu_1, mu_2, Sigma_1, Sigma_2)
    wasserstein_distance = norm(mu_1 - mu_2)^2 + trace(Sigma_1) + trace(Sigma_2) - 2 * trace((Sigma_1^0.5 * Sigma_2 * Sigma_1^0.5)^0.5);
end