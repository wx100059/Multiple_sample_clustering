plot(k,ave_MI_K_means);
hold on;
plot(k,ave_MI_spectral);
plot(k,ave_MI_KL);
% plot(k,ave_MI_KL_plus);
%plot(k,ave_MI_KL_spectral);
% plot(k,ave_MI_WA_spectral);
plot(k,ave_MI_BH);
hold off
% legend('k-means','spectral clustering','KL divergence', 'KL divergence++', 'Wasserstein distance','bhattacharya distance');
legend('k-means','spectral clustering','KL divergence', 'bhattacharya distance');