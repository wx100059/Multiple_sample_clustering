function assignment = random_assign(probability_distribution)
%randomly assign the point based on the given distribution 
%   �˴���ʾ��ϸ˵��
            % make it into a cumulative distribution
            culmulative_distribution = cumsum(probability_distribution);
            random_number = rand; % Your trials
            assignment = find(culmulative_distribution > random_number, 1, 'first');
end

