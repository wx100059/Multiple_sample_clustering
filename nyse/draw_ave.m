for i = 1:max_noise_level
    figure;
    plot(k, squeeze(ave_MI_K_means(:,i)));
    hold on;
    plot(k, squeeze(ave_MI_KL(:,i)));
    plot(k, squeeze(ave_MI_KL_plus(:,i)));
    xlabel('k');
    ylabel('normalized mutual information');
    hold off
    legend('k-means','KL divergence', 'KL divergence++');
    figure
    plot(k, squeeze(ave_MI_spectral(:,i)));
    hold on
    plot(k, squeeze(ave_MI_WA(:,i)));
    plot(k, squeeze(ave_MI_BH(:,i)));
    xlabel('k');
    ylabel('normalized mutual information');
    hold off
    legend('spectral clustering', 'Wasserstein distance','Bhattacharya distance');
end