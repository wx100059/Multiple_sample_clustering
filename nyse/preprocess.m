% this script does all the preprocess and Gaussian abstract distribution of
% each stocks from the origional data.
clear;
clc;
tic
% read data and build the index of different stocks
T = readtable('prices-split-adjusted.csv');
unique_T = unique(T(:,2));
T_index = cell(1,height(unique_T));
for i = 1:height(unique_T)
    T_index{i} = find(ismember(T(:,2), unique_T(i,1)));
end

% There are 501 stocks and we choose 500 of them to make sure all stocks
% has the same sample number
feature_begin = 3;
feature_end = 6;
d = feature_end - feature_begin +1;
total_estimated_Gaussian = struct('mean', zeros(1,5), 'covariance', zeros(d,d));
total_estimated_Gaussian(height(unique_T)) = struct('mean', zeros(1,5), 'covariance', zeros(d,d));
for i = 1:height(unique_T)
    total_estimated_Gaussian(i).mean = mean(log2(1 + table2array(T(T_index{i},feature_begin:feature_end))),1);
    total_estimated_Gaussian(i).covariance = cov(log2(1 + table2array(T(T_index{i},feature_begin:feature_end))));
end
save preprocess
toc