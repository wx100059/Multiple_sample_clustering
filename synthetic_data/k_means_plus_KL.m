function kl_assign = k_means_plus_KL(total_estimated_Gaussian, k, d, object_num, max_it)
%K_MEANS_PLUS_KL 此处显示有关此函数的摘要
%   此处显示详细说明
ave_mu = zeros(k,d);
ave_covariance = zeros(k,d,d);
initial_cluster_center = randi(object_num,1,k);
initial_cluster_center(1) = randi(object_num);
ave_mu(1,:) = total_estimated_Gaussian(initial_cluster_center(1)).mean;
ave_covariance(1,:,:) = total_estimated_Gaussian(initial_cluster_center(1)).covariance;
it_index = 0;
initial_assign = zeros(1,object_num);
it_count = 0;
record_assign = zeros(1,object_num);
        if k >= 2
            for j = 2:k
                temp_distance = zeros(1,j-1);       
                initial_distance = zeros(1, object_num);            
                for i = 1:object_num
                    % choose the smallest distance between the i-th point and all
                    % chosen cluster centers.
                    for p = 1:(j-1) % p is an indicator to show the nearest distance between one point to the chosen cluster center. 
                        if ismember(i,initial_cluster_center) == false
                            temp_distance(p) = kl_compute(total_estimated_Gaussian(i).mean, total_estimated_Gaussian(initial_cluster_center(p)).mean, total_estimated_Gaussian(i).covariance, total_estimated_Gaussian(initial_cluster_center(p)).covariance);
                           initial_distance(i) = min(temp_distance); 
                        end
                    end            
                end
                % randomly choose the initial cluster center in based on their
                % distance between the point and the nearest cluster centers
                probability_distribution = initial_distance.^2 / sum(initial_distance.^2);
                initial_cluster_center(j) = random_assign(probability_distribution);
                ave_mu(j,:) = total_estimated_Gaussian(initial_cluster_center(j)).mean;
                ave_covariance(j,:,:) = total_estimated_Gaussian(initial_cluster_center(j)).covariance;
            end
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

