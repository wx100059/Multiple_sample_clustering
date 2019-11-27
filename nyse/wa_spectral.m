function wa_spectral_assign = wa_spectral(total_estimated_Gaussian, k, object_num)
%UNTITLED 此处显示有关此函数的摘要
%   build the distance graph between distributions
       wa_graph = zeros(object_num, object_num);
       for i = 2:object_num
           for j = 1:(i-1)
               wa_graph(i,j) = wa_compute(total_estimated_Gaussian(i).mean, total_estimated_Gaussian(j).mean,total_estimated_Gaussian(i).covariance,total_estimated_Gaussian(j).covariance);
           end
       end
       
       for i = 1:object_num
           for j = (i+1):object_num
               wa_graph(i,j) = wa_graph(j,i);
           end
       end
% transform the distance graph into similarity graph
%     wa_graph_max = max(max(wa_graph));
%     normalized_wa_graph = 1 - wa_graph/wa_graph_max;
    % d = length(total_estimated_Gaussian(1).mean);
%     if d == 9
%         temp = 0;
%     end
%     if d ==10
%         record = 0;
%     end
    % sigma = 0.35 * d^2 - 0.4 * d -1;

    sigma = mean(mean(wa_graph));
    normalized_wa_graph = exp(-wa_graph.^2./(2*sigma^2));

       Z = spectral_clustering(normalized_wa_graph, k, 0);
       wa_spectral_assign = zeros(1,object_num);
       for i = 1:k
           index = find(Z(:,i));
           wa_spectral_assign(index) = i;
       end
end

