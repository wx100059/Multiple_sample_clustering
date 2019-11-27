function KL_spectral_assign = KL_spectral(total_estimated_Gaussian, k,object_num)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    X = zeros(object_num); % stores the KL divergence between each pair of distributions. 
    %sigma = 1; % the parameter of normalization
    for i = 1:object_num
        for j = 1:object_num
            X(i,j) = kl_compute(total_estimated_Gaussian(i).mean, total_estimated_Gaussian(j).mean, total_estimated_Gaussian(i).covariance, total_estimated_Gaussian(j).covariance);
        end
    end
    balanced_X = zeros(object_num); 
    for i = 1:object_num
        for j = 1:object_num
            balanced_X(i,j) = (X(i,j) + X(j,i))/2;
        end
    end
    d = length(total_estimated_Gaussian(1).mean);
%     sigma = 0.35 * d^2 - 0.4 * d -1;
    sigma = mean(mean(balanced_X));
    normalized_X = exp(-balanced_X.^2./(2*sigma^2));
% 	balanced_X_max = max(max(balanced_X));
%     normalized_X = 1 - balanced_X/balanced_X_max;
    KL_spectral = spectral_clustering(normalized_X, k ,0);
    KL_spectral_assign = zeros(1,object_num);
    for i = 1:k
        index = find(KL_spectral(:,i));
        KL_spectral_assign(index) = i;
    end
end

