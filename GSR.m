% This function is used to re-orthonormalize the column vectors of the matrix a.
% E - second-paradigm number of each column vector after orthogonalization
% b - re-orthonormalize result
function [b,E]=GSR(a)
    [m,n] = size(a);
    if(m<n)
        error('Row is smaller than column, cannot be calculated, please re-enter after transposing');
        return
    end
    b=zeros(m,n);
    b(:,1)=a(:,1);
    for i=2:n
        for j=1:i-1
            b(:,i)=b(:,i)-dot(a(:,i),b(:,j))/dot(b(:,j),b(:,j))*b(:,j);
        end
        b(:,i)=b(:,i)+a(:,i);
    end


    for k=1:n
        E(k)=norm(b(:,k));
        b(:,k)=b(:,k)/E(k);
    end
    E=E';
end