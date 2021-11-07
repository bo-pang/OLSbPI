
% This function vectorizes a symmetric matrix

function x = sm2vec(X)
[n,~] = size(X);
N = n*(n+1)/2;
x = zeros(N,1);
k = 1;
for i=1:n
    for j=i:n
        if i==j
            x(k) = X(i,j);
        else
            x(k) = sqrt(2)*X(i,j);
        end
        k = k+1;
    end
end
end