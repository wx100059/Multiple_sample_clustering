histogram(k_means_assign);
title('k means++');
figure
histogram(spectral_assign);
title('spectral_clustering');
figure
histogram(kl_assign);
title('KL divergence');
figure
histogram(kl_plus_assign);
title('KL divergence k means++');
figure
histogram(wa_spectral_assign);
title('Wasserstein distance');