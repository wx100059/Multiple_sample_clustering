function mean_vector = simplex_sample(d)
% generate the Gaussian mean vector from sampling the unit simplex
% the detailed provement and algorithm can be seen in 
% https://www.cs.cmu.edu/~nasmith/papers/smith+tromble.tr04.pdf
% d: the dimension of unit simplex
    M = intmax('int32');  
    x = randi([1, M], 1, d-1);
    x = sort(x);
    x = [0, x, M];
    mean_vector = zeros(1, d);
    for i = 1:d
        mean_vector(i) = x(i+1) - x(i);
    end
    mean_vector = mean_vector/double(M);
end

