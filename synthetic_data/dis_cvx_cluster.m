function [cvx_assign4] = dis_cvx_cluster(Gaussian,k)
object_num = length(Gaussian);
d  = length(Gaussian(1).mean);
lambda1 = 0.4;
lambda2 = 0.2;
sigma1 = 1;
sigma2 = 2;
w_mean = zeros(object_num,object_num);
w_cov = zeros(object_num,object_num);
temp = zeros(1,object_num);
I_mean = zeros(object_num, k);
I_cov = zeros(object_num , k);
for i = 1:object_num
    for j = 1:object_num
        temp(j) = exp(-(norm(Gaussian(i).mean - Gaussian(j).mean,2))/sigma1^2);
    end
    [B_mean, I_mean(i,:)] = maxk(temp, k);
%     w_mean(i,I_mean(i,:)) = B_mean;
    w_mean(i,I_mean(i,:)) =1;
end

for i = 1:object_num
    for j = 1:object_num
        temp(j) = exp(-norm((Gaussian(i).covariance(:) - Gaussian(j).covariance(:))/sigma2^2,2));
    end
    [B_cov, I_cov(i,:)] = maxk(temp, k);
%     w_cov(i,I_cov(i,:)) = B_cov;
    w_cov(i,I_cov(i,:)) = 1;
end


cvx_begin
    variable mean_vec(object_num, d)
    sum_fun1 = 0;

    
    for i = 1:object_num
        sum_fun1 = sum_fun1 + norm((mean_vec(i,:) - Gaussian(i).mean),2);
    end
    
    for i = 1:object_num
        for j = I_mean(i,:)
            sum_fun1 = sum_fun1 + lambda1 * w_mean(i,j) * norm(mean_vec(i,:) - mean_vec(j,:), 2);
        end
    end   
    minimize(sum_fun1)
cvx_end

cvx_begin
    variable mean_cov(object_num,d,d)
    sum_fun2 = 0;
    for i = 1:object_num
        sum_fun2 = sum_fun2 + norm(mean_cov(i,:)' - Gaussian(i).covariance(:), 2);
    end

    for i = 1:object_num
        for j = I_cov(i,:)
            sum_fun2 = sum_fun2 + lambda2 * w_cov(i,j) * norm(mean_cov(i,:) - mean_cov(j,:),2);
        end
    end
    minimize(sum_fun2)
cvx_end
estimated_Gaussian.mean = mean_vec(1,:);
estimated_Gaussian.covariance = squeeze(mean_cov(1,:,:));
total_estimated_Gaussian(1) = estimated_Gaussian;
total_estimated_Gaussian(object_num) = estimated_Gaussian;

for i = 1:object_num
    total_estimated_Gaussian(i).mean = mean_vec(i,:);
    total_estimated_Gaussian(i).covariance = squeeze(mean_cov(i,:,:));
end

% cvx_assign1 = wa_spectral(total_estimated_Gaussian, k, object_num);
% cvx_assign2 = bh_spectral(total_estimated_Gaussian, k, object_num);
% cvx_assign3 = k_means_KL(total_estimated_Gaussian, k, d, object_num, 100);
cvx_assign4 = k_means_plus_KL(total_estimated_Gaussian, k, d, object_num, 100);
dis_mat = zeros(object_num,object_num);
for i = 1:object_num
    for j = 1:i
        dis_mat(i,j) = wa_compute(mean_vec(i,:), mean_vec(j,:), squeeze(mean_cov(i,:,:)),squeeze(mean_cov(j,:,:)));
    end
end
for i = 1:(object_num-1)
    for j = (i+1):object_num
        dis_mat(i,j) = dis_mat(j,i);
    end
end


sigma = mean(mean(dis_mat));
normalized_dis_mat = exp(-dis_mat.^2./(2*sigma^2));
Z = spectral_clustering(normalized_dis_mat, k, 0);
cvx_spectral_assign = zeros(1,object_num);
 for i = 1:k
     index = find(Z(:,i));
     cvx_spectral_assign(index) = i;
 end