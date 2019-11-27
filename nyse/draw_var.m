for i = 1:max_noise_level
    figure;
    plot(k, squeeze(var_MI_K_means(:,i)));
    hold on;
    plot(k, squeeze(var_MI_KL(:,i)));
    plot(k, squeeze(var_MI_KL_plus(:,i)));
    hold off
    legend('k-means','KL divergence', 'KL divergence++');
    figure
    plot(k, squeeze(var_MI_spectral(:,i)));
    hold on
    plot(k, squeeze(var_MI_WA(:,i)));
    plot(k, squeeze(var_MI_BH(:,i)));
    hold off
    legend('spectral clustering', 'Wasserstein distance','Bhattacharya distance');
end