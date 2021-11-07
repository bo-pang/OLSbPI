
% This function implements the operation $svec(xx^T)$

function X = kronv(x)
len = length(x);
X = [];
for i=1:len
    for j=i:len
        if i==j
            X(end+1) = x(i)*x(j);
        else
            X(end+1) = sqrt(2)*x(i)*x(j);
        end
    end
end
X = X';
end