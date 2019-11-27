function spectral_assign = spectral(X, k, object_num)
%SPECTRAL 此处显示有关此函数的摘要
%   此处显示详细说明
        X_graph = zeros(object_num, object_num);
        for i = 1:object_num
            for j = 1:object_num
                X_graph(i,j) = norm(X(:,i) - X(:,j))^2;
            end
        end

        X_max = max(max(X_graph));
        for i = 1:object_num
            for j = 1:object_num
                X_graph(i,j) = 1 - X_graph(i,j)/X_max;
            end
        end
        
       X_spectral = spectral_clustering(X_graph, k, 0);
       spectral_assign = zeros(1,object_num);
       for i = 1:k
           index = find(X_spectral(:,i));
           spectral_assign(index) = i;
       end
end

