clear;
clc;
tic
k = 5; % k represents the number of clusters
sample_num = 30; %sample number from each distribution
object_num = 200; % obejct number. Each obejct is a set of sample from the same Gaussian distribution
real_cluster = zeros(1, object_num); % the ground truth used to compare the clustering accuracy.
loop_time = 1000;
max_it = 30;
%noise_level = 1.41;
d_interval = 7; % d_interval represents the range of d
d_begin = 4; % where d starts

%store the mutual information between assignment and ground truth
MI_K_means = zeros(d_interval, loop_time);
MI_spectral = zeros(d_interval, loop_time);
MI_KL = zeros(d_interval, loop_time);
MI_KL_plus = zeros(d_interval, loop_time);
%MI_KL_spectral = zeros(d_interval, loop_time);
MI_WA_spectral = zeros(d_interval, loop_time);
MI_BH = zeros(d_interval, loop_time);

%compute the average mutual information after 1000 loop times
ave_MI_K_means = zeros(1, d_interval);
ave_MI_spectral = zeros(1, d_interval);
ave_MI_KL = zeros(1, d_interval);
ave_MI_KL_plus = zeros(1, d_interval);
%ave_MI_KL_spectral = zeros(1, d_interval);
ave_MI_WA_spectral = zeros(1, d_interval);
ave_MI_BH = zeros(1, d_interval);

var_MI_K_means = zeros(1, d_interval);
var_MI_spectral = zeros(1, d_interval);
var_MI_KL = zeros(1, d_interval);
var_MI_KL_plus = zeros(1, d_interval);
% var_MI_KL_spectral = zeros(1, k_interval);
var_MI_WA_spectral = zeros(1, d_interval);
var_MI_BH = zeros(1, d_interval);

%create object_num number of estimated distributions from k predefined
%distirbutions
for d = d_begin:(d_begin+d_interval-1)
    total_sample = zeros(object_num, sample_num, d); % total_sample represent the a sample set that each obejct represents 30 sample from one Gaussian distribution and 200 obejcts is considered.
    object = struct('mean', ones(1,d), 'covariance', diag(ones(1,d)));% create a structure that represent one Gaussian distribution.
    total_Gaussian = object;% create a vecotr that contains k Gaussian distribution.
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
        % make sample from Gaussian distribution
        d
        loop_index
        for i = 1:object_num
            total_sample(i, : ,:) = mvnrnd(total_Gaussian(mod(i-1, k) + 1).mean, total_Gaussian(mod(i-1, k) + 1).covariance, sample_num);
            real_cluster(i) = mod(i-1, k) + 1;
        end

        estimated_Gaussian = struct('mean',zeros(1:d), 'covariance', zeros(d)); % store the Gasusain reconstructed from the sample
        total_estimated_Gaussian = estimated_Gaussian;  % store the all the Gaussian reconstructed from the sample
        total_estimated_Gaussian(object_num)= estimated_Gaussian;
        % reconstruct the Gaussian distribution from the sample
        % This is an unbiased estimation
        for i = 1:object_num
            total_estimated_Gaussian(i).mean = mean(squeeze(total_sample(i,:,:)),1);
            total_estimated_Gaussian(i).covariance = cov(squeeze(total_sample(i,:,:)));
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

%         % use symmetric KL divergence matrix to do spectral clustering
%         kl_spectral_assign = KL_spectral(total_estimated_Gaussian, k, object_num);
        
       % use wasserstein distance based spectral clustering
        wa_spectral_assign = wa_spectral(total_estimated_Gaussian, k, object_num);
        
       % use bhattacharya distance based spectral clustering
        bh_spectral_assign = bh_spectral(total_estimated_Gaussian, k, object_num);

        MI_K_means(d-d_begin+1,loop_index) = nmi(real_cluster, k_means_assign);
        MI_spectral(d-d_begin+1, loop_index) = nmi(real_cluster, spectral_assign);
        MI_KL(d-d_begin+1, loop_index) = nmi(real_cluster, kl_assign);
%         MI_KL_spectral(d-d_begin+1, loop_index) = nmi(real_cluster, kl_spectral_assign);
        MI_KL_plus(d-d_begin+1, loop_index) = nmi(real_cluster, kl_plus_assign);
        MI_WA_spectral(d-d_begin+1, loop_index) = nmi(real_cluster, wa_spectral_assign);
        MI_BH(d-d_begin+1, loop_index) = nmi(real_cluster, bh_spectral_assign);        

        
    end
    ave_MI_KL(d-d_begin+1) = mean(MI_KL(d-d_begin+1, :));
%     ave_MI_KL_spectral(d-d_begin+1) = mean(MI_KL_spectral(d-d_begin+1, :));
    ave_MI_KL_plus(d-d_begin+1) = mean(MI_KL_plus(d-d_begin+1, :));
    ave_MI_WA_spectral(d-d_begin+1) = mean(MI_WA_spectral(d-d_begin+1, :));
    ave_MI_K_means(d-d_begin+1) = mean(MI_K_means(d-d_begin+1, :));
    ave_MI_spectral(d-d_begin+1) = mean(MI_spectral(d-d_begin+1, :));
    ave_MI_BH(d-d_begin+1) = mean(MI_BH(d-d_begin+1, :));
    %ave_MI_WA_k_means_plus(d-d_begin+1) = mean(MI_WA_k_means_plus(d-d_begin+1, :));
    
    var_MI_K_means(d-d_begin+1) = var(MI_K_means(d-d_begin+1, :));
    var_MI_spectral(d-d_begin+1) = var(MI_spectral(d-d_begin+1, :));
%     var_MI_KL_spectral(k-k_begin+1) = var(MI_KL_spectral(d-d_begin+1, :));
    var_MI_KL(d-d_begin+1) = var(MI_KL(d-d_begin+1, :));
    var_MI_KL_plus(d-d_begin+1) = var(MI_KL_plus(d-d_begin+1, :));
    var_MI_WA_spectral(d-d_begin+1) = var(MI_WA_spectral(d-d_begin+1, :));
    var_MI_BH(d-d_begin+1) = var(MI_BH(d-d_begin+1, :));
end
run_time = toc;
% save 1000_k_5.mat;
toc