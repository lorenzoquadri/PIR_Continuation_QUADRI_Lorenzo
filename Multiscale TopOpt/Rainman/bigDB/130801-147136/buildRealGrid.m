%parpool('local',24)

fileID=fopen('gridNoNoise.txt','w');
gridDim=[4 4 4 4 4 4];
gridSize=gridDim(1)*gridDim(2)*gridDim(3)*gridDim(4)*gridDim(5)*gridDim(6);
bigGrid=zeros(gridSize,6);
for i=1:gridDim(1)
    for j=1:gridDim(2)
        for k=1:gridDim(3)
            for l=1:gridDim(4)
                for m=1:gridDim(5)
                    for n=1:gridDim(6)
                        x1=i/(gridDim(1)+1); %0 and 1 values are excluded
                        x2=(j-1)/(gridDim(2)-1); %0 and 1 values are included
                        %x2=min(abs(x2),2-x2); %put value in [0,1]
                        x3=(k-1)/(gridDim(3)-1); %0 and 1 values are included
                        %x3=min(abs(x3),2-x3); %put value in [0,1]
                        x4=(l-1)/(gridDim(4)-1); %0 and 1 values are included
                        x5=(m-1)/(gridDim(5)-1); %0 and 1 values are included
                        x6=(n-1)/(gridDim(6)-1); %0 and 1 values are included
                        % bigGrid((i-1)*gridDim(2)*gridDim(3)+(j-1)*gridDim(3)+k,:)=[x1 x2 x3];
                        bigGrid((i-1)*gridDim(2)*gridDim(3)*gridDim(4)*gridDim(5)*gridDim(6)+(j-1)*gridDim(3)*gridDim(4)*gridDim(5)*gridDim(6)+(k-1)*gridDim(4)*gridDim(5)*gridDim(6)+(l-1)*gridDim(5)*gridDim(6)+(m-1)*gridDim(6)+n,:)=[x1 x2 x3 x4 x5 x6];
                    end
                end
            end
        end
    end
end

for i=1:size(bigGrid,1)
    fprintf(fileID,'%14.10f',bigGrid(i,:));
    fprintf(fileID,'\n');
end
fclose('all')
