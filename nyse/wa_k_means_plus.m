function wa_assign = wa_k_means_plus(total_estimated_Gaussian, k, d, object_num, max_it)
%WA_K_MEANS_PLUS 此处显示有关此函数的摘要
%   此处显示详细说明
    %     % use k-means++ to assign the initial Wasserstein distance clusters
    %     % randomly choose the first cluster center
initial_cluster_center = zeros(1,k);
ave_mu = zeros(k,d);
ave_covariance = zeros(k,d,d);
initial_assign = zeros(1,object_num);
record_assign = zeros(1,object_num);
        initial_cluster_center(1) = randi(object_num);
        ave_mu(1,:) = total_estimated_Gaussian(initial_cluster_center(1)).mean;
        ave_covariance(1,:,:) = total_estimated_Gaussian(initial_cluster_center(1)).covariance;
        if k >= 2
            for j = 2:k
                temp_distance = zeros(1,j-1);       
                initial_distance = zeros(1, object_num);            
                for i = 1:object_num
                    % choose the smallest distance between the i-th point and all
                    % chosen cluster centers.
                    for p = 1:(j-1) % p is an indicator to show the nearest distance between one point to the chosen cluster center. 
                        if ismember(i,initial_cluster_center) == false
                            temp_distance(p) =wa_compute(total_estimated_Gaussian(i).mean, total_estimated_Gaussian(initial_cluster_center(p)).mean, total_estimated_Gaussian(i).covariance, total_estimated_Gaussian(initial_cluster_center(p)).covariance);
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
        temp_wa_distance = zeros(1,k);
         for i = 1:object_num
            for j = 1:k
                temp_wa_distance(j) = wa_compute(total_estimated_Gaussian(i).mean, ave_mu(j,:), total_estimated_Gaussian(i).covariance, squeeze(ave_covariance(j,:,:)));
            end
            index = find(temp_wa_distance == min(temp_wa_distance));
            initial_assign(i) = index(1);
         end
         
        % use wasserstein distance based algorithm to cluster the algorithm
        wa_assign = initial_assign;
        it_count = 0;
        while (isequal(record_assign, wa_assign) == false) && (it_count < max_it)
            % compute the mean vector of the Gaussian barycenter 
            it_count = it_count + 1;
            for j = 1:k
                index_j = find(wa_assign == j);
                temp_mu = zeros(length(index_j),d);
                if(isempty(index_j))
                    continue;
                end
                for i = 1:length(index_j)
                    temp_mu(i,:) = total_estimated_Gaussian(index_j(i)).mean;
                end
                ave_mu(j,:) = mean(temp_mu, 1);
            end
    
            % compute the covariance matrix of the Gaussian barycenter
            for j = 1:k
                index_j = find(wa_assign == j);
                if(isempty(index_j))
                    continue;
                end    
                ave_covariance(j,:,:) = wasserstein_barycenter(total_estimated_Gaussian, index_j, d);
            end
            % compute Wasserstein distance between estimated Gaussian distribution and
            % the clustered average Gaussian distribution
            record_assign = wa_assign;
            for i = 1:object_num
                for j = 1:k
                    temp_wa_distance(j) = wa_compute(total_estimated_Gaussian(i).mean, ave_mu(j,:), total_estimated_Gaussian(i).covariance, squeeze(ave_covariance(j,:,:)));
                end
                index = find(temp_wa_distance == min(temp_wa_distance));
                wa_assign(i) = index(1);
            end    
        end
end

