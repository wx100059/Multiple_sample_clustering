clear;
clc;
tic
d = 4; % d represents the dimension of Gaussian multivariate
sample_num = 30; %sample number from each distribution
object_num = 50; % obejct number. Each obejct is a set of sample from the same Gaussian distribution
total_sample = zeros(object_num, sample_num, d); % total_sample represent the a sample set that each obejct represents 30 sample from one Gaussian distribution and 200 obejcts is considered.
object = struct('mean', ones(1,d), 'covariance', diag(ones(1,d)));% create a structure that represent one Gaussian distribution.
total_Gaussian = object;% create a vecotr that contains k Gaussian distribution.
k_interval = 1; % k_interval represents the range of k
k_begin = 5; % where k starts
real_cluster = zeros(k_interval, object_num); % the ground truth used to compare the clustering accuracy.
loop_time = 10;
%noise_level = 1;
max_it = 30;
MI_K_means = zeros(k_interval, loop_time);
MI_spectral = zeros(k_interval, loop_time);
MI_KL = zeros(k_interval, loop_time);
MI_KL_plus = zeros(k_interval, loop_time);
%MI_KL_spectral = zeros(k_interval, loop_time);
MI_WA_spectral = zeros(k_interval, loop_time);
MI_BH = zeros(k_interval, loop_time);
MI_CVX = zeros(k_interval, loop_time);

ave_MI_K_means = zeros(1, k_interval);
ave_MI_spectral = zeros(1, k_interval);
ave_MI_KL = zeros(1, k_interval);
ave_MI_KL_plus = zeros(1, k_interval);
% ave_MI_KL_spectral = zeros(1, k_interval);
ave_MI_WA_spectral = zeros(1, k_interval);
ave_MI_BH = zeros(1, k_interval);
ave_MI_CVX = zeros(1, k_interval);

var_MI_K_means = zeros(1, k_interval);
var_MI_spectral = zeros(1, k_interval);
var_MI_KL = zeros(1, k_interval);
var_MI_KL_plus = zeros(1, k_interval);
% var_MI_KL_spectral = zeros(1, k_interval);
var_MI_WA_spectral = zeros(1, k_interval);
var_MI_BH = zeros(1, k_interval);
var_MI_CVX = zeros(1, k_interval);

for k = k_begin:(k_begin + k_interval - 1)
    if k > 1  % assign value for the Gaussian distribution, the mean vector is randomly sampled from
        % unit simplex and the covariance matrix from the set of matrices with
        % eigenvalues 1,2,...,d.
        total_Gaussian(k)=object;
        for iterator = 1:k
           total_Gaussian(iterator)=object;
           total_Gaussian(iterator).mean = simplex_sample(d);
           temp = RandOrthMat(d);
           total_Gaussian(iterator).covariance = temp * diag(1:d) * temp';
        end
    end

    for loop_index = 1:loop_time
        k
        loop_index
        % make sample from Gaussian distribution
        for i = 1:object_num
            total_sample(i, : ,:) = mvnrnd(total_Gaussian(mod(i-1, k) + 1).mean, total_Gaussian(mod(i-1, k) + 1).covariance, sample_num);
            real_cluster(k - k_begin + 1,i) = mod(i-1, k) + 1;
        end

        estimated_Gaussian = struct('mean',zeros(1:d), 'covariance', zeros(d)); % store the Gasusain reconstructed from the sample
%         total_estimated_Gaussian = estimated_Gaussian;  % store the all the Gaussian reconstructed from the sample
%         total_estimated_Gaussian(object_num)= estimated_Gaussian;
        % reconstruct the Gaussian distribution from the sample
        % This is an unbiased estimation
        for i = 1:object_num
            total_estimated_Gaussian(k-k_begin+1,i).mean = mean(squeeze(total_sample(i,:,:)),1);
            total_estimated_Gaussian(k-k_begin+1,i).covariance = cov(squeeze(total_sample(i,:,:)));
        end

        % abstract the first order information
        X = zeros(d,object_num);
        for i = 1:object_num
            X(:,i) = total_estimated_Gaussian(k-k_begin+1,i).mean;
        end
        
        % use k-means++ to do the clustering based on the first order information
        k_means_assign = (kmeans(X', k))';

        % Do spectral clustering use the first order information
        spectral_assign = spectral(X, k, object_num);

        % use k-means based KL divergence clusters 
        kl_assign = k_means_KL(total_estimated_Gaussian(k-k_begin+1,:), k, d, object_num, max_it);
        
        % use k-means++ based KL divergence clusters
        kl_plus_assign = k_means_plus_KL(total_estimated_Gaussian(k-k_begin+1,:), k, d, object_num, max_it);

%         % use convex cluster 
        cvx_cluster_assign = dis_cvx_cluster(total_estimated_Gaussian,k);
%         % use symmetric KL divergence matrix to do spectral clustering
%         kl_spectral_assign = KL_spectral(total_estimated_Gaussian, k, object_num);
        
       % use wasserstein distance based sepectral clustering
        wa_spectral_assign = wa_spectral(total_estimated_Gaussian(k-k_begin+1,:), k, object_num);
         
       % use bhattacharya distance based spectral clustering
        bh_spectral_assign = bh_spectral(total_estimated_Gaussian(k-k_begin+1,:), k, object_num);
        
        MI_K_means(k-k_begin+1,loop_index) = nmi(real_cluster(k - k_begin + 1, :), k_means_assign);
        MI_spectral(k-k_begin+1, loop_index) = nmi(real_cluster(k - k_begin + 1, :), spectral_assign);
        MI_KL(k-k_begin+1, loop_index) = nmi(real_cluster(k - k_begin + 1, :), kl_assign);
%         MI_KL_spectral(k-k_begin+1, loop_index) = nmi(real_cluster, kl_spectral_assign);
        MI_KL_plus(k-k_begin+1, loop_index) = nmi(real_cluster(k - k_begin + 1, :), kl_plus_assign);
        MI_WA_spectral(k-k_begin+1, loop_index) = nmi(real_cluster(k - k_begin + 1, :), wa_spectral_assign);
        MI_BH(k-k_begin+1, loop_index) = nmi(real_cluster(k - k_begin + 1, :), bh_spectral_assign);  
        MI_CVX(k-k_begin+1, loop_index) = nmi(real_cluster(k - k_begin + 1, :), cvx_cluster_assign);  
    end

    ave_MI_K_means(k-k_begin+1) = mean(MI_K_means(k-k_begin+1, :));
    ave_MI_spectral(k-k_begin+1) = mean(MI_spectral(k-k_begin+1, :));
%     ave_MI_KL_spectral(k-k_begin+1) = mean(MI_KL_spectral(k-k_begin+1, :));
    ave_MI_KL(k-k_begin+1) = mean(MI_KL(k-k_begin+1, :));
    ave_MI_KL_plus(k-k_begin+1) = mean(MI_KL_plus(k-k_begin+1, :));
    ave_MI_WA_spectral(k-k_begin+1) = mean(MI_WA_spectral(k-k_begin+1, :));
    ave_MI_BH(k-k_begin+1) = mean(MI_BH(k-k_begin+1, :));
    ave_MI_CVX(k-k_begin+1) = mean(MI_CVX(k-k_begin+1, :));
    
    var_MI_K_means(k-k_begin+1) = var(MI_K_means(k-k_begin+1, :));
    var_MI_spectral(k-k_begin+1) = var(MI_spectral(k-k_begin+1, :));
%     var_MI_KL_spectral(k-k_begin+1) = var(MI_KL_spectral(k-k_begin+1, :));
    var_MI_KL(k-k_begin+1) = var(MI_KL(k-k_begin+1, :));
    var_MI_KL_plus(k-k_begin+1) = var(MI_KL_plus(k-k_begin+1, :));
    var_MI_WA_spectral(k-k_begin+1) = var(MI_WA_spectral(k-k_begin+1, :));
    var_MI_BH(k-k_begin+1) = var(MI_BH(k-k_begin+1, :));
    var_MI_CVX(k-k_begin+1) = var(MI_CVX(k-k_begin+1, :));
end
run_time = toc;
%save 1000_d_7.mat;
toc