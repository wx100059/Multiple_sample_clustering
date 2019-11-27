%compute the kl divergence between two Gaussian distributions
function D_kl = kl_compute(mu_1, mu_2, Sigma_1, Sigma_2)
dim = size(Sigma_1);
D_kl = 1 / 2 * (log2(det(Sigma_2) / det(Sigma_1)) - dim(1) + trace(Sigma_2^(-1) * Sigma_1) + (mu_2 - mu_1) * Sigma_2^(-1) * (mu_2 - mu_1)');
end