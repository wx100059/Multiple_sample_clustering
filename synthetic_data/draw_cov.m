plot(k,var_MI_K_means);
hold on;
plot(k,var_MI_spectral);
plot(k,var_MI_KL);
plot(k,var_MI_KL_plus);
%plot(k,ave_MI_KL_spectral);
plot(k,var_MI_WA_spectral);
plot(k,var_MI_BH);
hold off
legend('k-means','spectral clustering','KL divergence', 'KL divergence++', 'Wasserstein distance','bhattacharya distance');
%legend('KL divergence', 'KL divergence++','KL spectral', 'Wasserstein distance','bhattacharya distance');
xlabel('d');
ylabel('variance');