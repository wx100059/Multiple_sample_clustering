function bh_spectral_assign = bh_spectral(total_estimated_Gaussian, k, object_num)
%UNTITLED 此处显示有关此函数的摘要
%   build the bhattacharya distance graph between distributions
       bh_graph = zeros(object_num, object_num);
       for i = 2:object_num
           for j = 1:(i-1)
               bh_graph(i,j) = bhattacharya_compute(total_estimated_Gaussian(i).mean, total_estimated_Gaussian(j).mean,total_estimated_Gaussian(i).covariance,total_estimated_Gaussian(j).covariance);
           end
       end
% d = length(total_estimated_Gaussian(1).mean);
% transform the distance graph into similarity graph
%     bh_graph_max = max(max(bh_graph));
%     normalized_bh_graph  = 1 - bh_graph/bh_graph_max;
%     sigma = 0.25 * d^2 - 0.5 * d -1.9 + 0.1 * k;
       for i = 1:object_num
           for j = (i+1):object_num
               bh_graph(i,j) = bh_graph(j,i);
           end
       end
    sigma = mean(mean(bh_graph));
    normalized_bh_graph = exp(-bh_graph.^2./(2*sigma^2));
%     if k == 8
%         temp = 0;
%     end
%     if d == 9
%         temp = 0;
%     end
%     if d ==10
%         record = 0;
%     end
       Z = spectral_clustering(normalized_bh_graph, k, 0);
       bh_spectral_assign = zeros(1,object_num);
       for i = 1:k
           index = find(Z(:,i));
           bh_spectral_assign(index) = i;
       end
end


