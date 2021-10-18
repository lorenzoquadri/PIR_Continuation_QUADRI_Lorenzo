function [matrix]=array2matrix(array)
dim=(-1+sqrt(1+4*2*length(array)))/2;
matrix=zeros(dim,dim);

loop=1;
for i=1:size(matrix,1)
    for j=1:size(matrix,2)
        if j>=i
            matrix(i,j)=array(loop);
            loop=loop+1;
        end
    end
end

% making the matrix symmetric
for i=1:size(matrix,1)
    for j=1:size(matrix,2)
        if j>=i
            matrix(j,i)=matrix(i,j);
        end
    end
end