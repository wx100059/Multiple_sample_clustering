% utilize KL-divergence, KL-divergence k means++, Wasserstein distance to
% do the clustering of stocks.
clear;
clc;
tic

load preprocess.mat
object_num = height(unique_T); % obejct number. Each obejct is a set of sample from the same Gaussian distribution
loop_time = 100;
max_it = 30; % biggest iteration time
max_noise_level = 3;
k_interval = 8 ; % k_interval represents the range of k
k_begin = 2; % where k starts

MI_K_means = zeros(k_interval, loop_time, max_noise_level);
MI_spectral = zeros(k_interval, loop_time, max_noise_level);
MI_KL = zeros(k_interval, loop_time, max_noise_level);
MI_KL_plus = zeros(k_interval, loop_time, max_noise_level);
MI_WA_spectral = zeros(k_interval, loop_time, max_noise_level);
MI_BH = zeros(k_interval, loop_time, max_noise_level);

ave_MI_KL_plus = zeros(k_interval, max_noise_level);
ave_MI_WA = zeros(k_interval, max_noise_level);
ave_MI_K_means = zeros(k_interval, max_noise_level);
ave_MI_spectral = zeros(k_interval, max_noise_level);
ave_MI_KL = zeros(k_interval, max_noise_level);
ave_MI_BH = zeros(k_interval, max_noise_level);

var_MI_KL_plus = zeros(k_interval, max_noise_level);
var_MI_WA = zeros(k_interval, max_noise_level);
var_MI_K_means = zeros(k_interval, max_noise_level);
var_MI_spectral = zeros(k_interval, max_noise_level);
var_MI_KL = zeros(k_interval, max_noise_level);
var_MI_BH = zeros(k_interval, max_noise_level);

for k = k_begin:(k_begin+k_interval-1)
    for loop_index = 1:loop_time
        % compute the ground truth
        initial_cluster_center = zeros(1,k); % choose k samples as the initial cluster center
        initial_assign = zeros(1,object_num);
        record_assign = zeros(1,object_num); % watch if the k-means algorithm converges
        ave_mu = zeros(k,d); % row of ave_mu stores the mean vector of on Gaussian distribution
        ave_covariance = zeros(k,d,d); % the average covariance matrix of k Gaussian distribution clusters.        
        X = zeros(d,object_num);
        for i = 1:object_num
           X(:,i) = total_estimated_Gaussian(i).mean;
        end

        % use k-means++ to do the clustering based on the first order information
        k_means_assign_ground_truth = (kmeans(X', k))';

        % Do spectral clustering use the first order information
        spectral_assign_ground_truth = spectral(X, k, object_num);

        % use k-means based KL divergence clusters 
        kl_assign_ground_truth = k_means_KL(total_estimated_Gaussian, k, d, object_num, max_it);

        % use k-means++ based KL divergence clusters
        kl_plus_assign_ground_truth = k_means_plus_KL(total_estimated_Gaussian, k, d, object_num, max_it);

       % use wasserstein distance based sepectral clustering
        wa_spectral_assign_ground_truth = wa_spectral(total_estimated_Gaussian, k, object_num);
   
       % use wasserstein distance based sepectral clustering
        bh_spectral_assign_ground_truth = bh_spectral(total_estimated_Gaussian, k, object_num);
           
        for noise_level = 1:max_noise_level
            for i = 1:height(unique_T)
                Gaussian_matrix = normrnd(0,noise_level,length(T_index{i}),d);
                total_estimated_Gaussian(i).mean = real(mean(log2(table2array(T(T_index{i},feature_begin:feature_end))+Gaussian_matrix),1));
                total_estimated_Gaussian(i).covariance = real(cov(log2(table2array(T(T_index{i},feature_begin:feature_end))+Gaussian_matrix)));
            end
            
            % abstract the first order information
            
            X = zeros(d,object_num);
            for i = 1:object_num
                X(:,i) = total_estimated_Gaussian(i).mean;
            end
            
            % use k-means++ to do the clustering based on the first order information
            k_means_assign = (kmeans(X', k))';

            % Do spectral clustering use the first order information
            spectral_assign = spectral(X, k, object_num);

            % use k-means based KL divergence clusters 
            kl_assign = k_means_KL(total_estimated_Gaussian, k, d, object_num, max_it);

            % use k-means++ based KL divergence clusters
            kl_plus_assign = k_means_plus_KL(total_estimated_Gaussian, k, d, object_num, max_it);

           % use wasserstein distance based sepectral clustering
            wa_spectral_assign = wa_spectral(total_estimated_Gaussian, k, object_num);

           % use bhattacharya distance based spectral clustering
            bh_spectral_assign = bh_spectral(total_estimated_Gaussian, k, object_num);
             
            MI_K_means(k-k_begin+1,loop_index, noise_level) = nmi(k_means_assign_ground_truth, k_means_assign);
            MI_spectral(k-k_begin+1, loop_index, noise_level) = nmi(spectral_assign_ground_truth, spectral_assign);
            MI_KL(k-k_begin+1, loop_index, noise_level) = nmi(kl_assign_ground_truth, kl_assign);
            MI_KL_plus(k-k_begin+1, loop_index, noise_level) = nmi(kl_plus_assign_ground_truth, kl_plus_assign);
            MI_WA_spectral(k-k_begin+1, loop_index, noise_level) = nmi(wa_spectral_assign_ground_truth, wa_spectral_assign);
            MI_BH(k-k_begin+1, loop_index, noise_level) = nmi(bh_spectral_assign_ground_truth, bh_spectral_assign);
                       
        end
    end

    ave_MI_KL = squeeze(mean(MI_KL, 2));
    ave_MI_KL_plus = squeeze(mean(MI_KL_plus, 2));
    ave_MI_WA = squeeze(mean(MI_WA_spectral, 2));
    ave_MI_K_means = squeeze(mean(MI_K_means, 2));
    ave_MI_spectral = squeeze(mean(MI_spectral, 2));
    ave_MI_BH = squeeze(mean(MI_BH, 2));
    
    var_MI_K_means = squeeze(var(MI_K_means, 0 ,2));
    var_MI_spectral = squeeze(var(MI_spectral, 0, 2));
    var_MI_KL = squeeze(var(MI_KL, 0, 2));
    var_MI_KL_plus = squeeze(var(MI_KL_plus, 0, 2));
    var_MI_WA = squeeze(var(MI_WA_spectral, 0, 2));
    var_MI_BH = squeeze(var(MI_BH, 0, 2));
end
k = k_begin:(k_begin+k_interval-1);
% save price_split_adjust;
run_time = toc;
toc