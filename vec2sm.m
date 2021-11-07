
% This function transfers a vector to the corresponding symmetric matrix  

function X = vec2sm(x,n)
X = zeros(n);
num = flip(1:n);
for i=1:n
    index = 0;
    for k=1:i-1
        index = index+num(k);
    end
    for j=0:n-i
        if j~=0
            X(i,i+j)=x(index+j+1)/sqrt(2);
            X(j+i,i)=X(i,j+i);
        else
            X(i,j+i)=x(index+j+1);
        end
    end
end
end