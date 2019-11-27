% compute the wasserstein barycenter of a set of Gaussian distributions
function ave_covariance = wasserstein_barycenter(total_estimated_Gaussian, index_j, d)
    % build the initial covariance matrix
    ave_covariance = diag(ones(1,d));
    record_covariance = zeros(d);
    counter = 0;
    % use iterate algorithm to aproximate the covariance matrix of wasserstein barycenter
    while sum(sum(record_covariance - ave_covariance).^2) > 0.5
        temp_covariance = zeros(d);
        for i = 1:length(index_j)
            temp_covariance = (temp_covariance + ave_covariance^0.5 * total_estimated_Gaussian(index_j(i)).covariance * ave_covariance^0.5)^0.5;
        end
        temp_covariance = 1 / length(index_j) * temp_covariance;
        record_covariance = ave_covariance;
        ave_covariance = ave_covariance^(-0.5) * temp_covariance^2 * ave_covariance^(-0.5);
        counter = counter + 1;
    end
end