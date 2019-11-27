function bh_distance = bhattacharya_compute(mu_1, mu_2, Sigma_1, Sigma_2)
%BHATTACHARYA_DISTANCE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
Sigma = (Sigma_1 + Sigma_2)/2;
bh_distance = 1/8*(mu_1 - mu_2) * Sigma^(-1) * (mu_1 - mu_2)' + 0.5 * log(det(Sigma) / sqrt(det(Sigma_1) * det(Sigma_2)));
end

