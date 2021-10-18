% Transform a symmetric matrix into a vector reading it line by line
function [array]=matrix2array(matrix)

if size(matrix,1)~=size(matrix,2)
    error("this is not a square matrix")
end

dim=sum(size(matrix,1):1);
array=zeros(1,dim);

loop=1;
for i=1:size(matrix,1)
    for j=1:size(matrix,2)
        if j>=i
            array(loop)=matrix(i,j);
            loop=loop+1;
        end
    end
end


end