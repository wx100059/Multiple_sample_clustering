function kl_assign = k_means_KL(total_estimated_Gaussian, k, d, object_num, max_it)
%K_MEANS_KL 此处显示有关此函数的摘要
%   此处显示详细说明
initial_cluster_center = randi(object_num,1,k);
ave_mu = zeros(k,d);
ave_covariance = zeros(k,d,d);
initial_assign = zeros(1,object_num);
record_assign = zeros(1,object_num);
it_index = 0;
for j = 1:k
    ave_mu(j,:) = total_estimated_Gaussian(initial_cluster_center(j)).mean;
    ave_covariance(j,:,:) = total_estimated_Gaussian(initial_cluster_center(j)).covariance;
end

        % build the initial cluster assignment
        temp_kl_divergence = zeros(1,k);
         for i = 1:object_num
            for j = 1:k
                temp_kl_divergence(j) = kl_compute(total_estimated_Gaussian(i).mean, ave_mu(j,:), total_estimated_Gaussian(i).covariance, squeeze(ave_covariance(j,:,:)));
            end
            index = find(temp_kl_divergence == min(temp_kl_divergence));
            initial_assign(i) = index(1);
         end

         % use kl divergence based cluster algorithm to cluster the gaussian
         % distribution.
        kl_assign = initial_assign;
       while (isequal(record_assign, kl_assign) == false) && (it_index < max_it)
            % update the average mean vecotr for each Gaussian cluster
            it_index = it_index + 1;
            for j = 1:k
                index_j = find(kl_assign == j);
                temp_mu = zeros(length(index_j),d);
                for i = 1:length(index_j)
                    temp_mu(i,:) = total_estimated_Gaussian(index_j(i)).mean;
                end
                ave_mu(j,:) = mean(temp_mu, 1);
            end

            % update the covariance matrix
            for j = 1:k
                index_j = find(kl_assign == j);
                temp_covariance = zeros(length(index_j),d,d);
                for i = 1:length(index_j)
                    temp_covariance(i,:,:) = total_estimated_Gaussian(index_j(i)).covariance + (total_estimated_Gaussian(index_j(i)).mean - ave_mu(j))' * (total_estimated_Gaussian(index_j(i)).mean - ave_mu(j));
                end
                ave_covariance(j,:,:) = squeeze(mean(temp_covariance,1));
            end
            % compute KL-divergence between estimated Gaussian distribution and
            % cluster average Gaussian distribution
            record_assign = kl_assign;
            for i = 1:object_num
                for j = 1:k
                    temp_kl_divergence(j) = kl_compute(total_estimated_Gaussian(i).mean, ave_mu(j,:), total_estimated_Gaussian(i).covariance, squeeze(ave_covariance(j,:,:)));
                end
                index = find(temp_kl_divergence == min(temp_kl_divergence));
                kl_assign(i) = index(1);
            end
      end
end

