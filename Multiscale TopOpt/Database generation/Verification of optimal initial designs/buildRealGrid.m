clear
clc
close all

density=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
angle1=[0 1/3 2/3 1];
angle2=[0 1/3 2/3 1];
angle3=[0 1/3 2/3 1];
cubicity21=[0 1/3 2/3 1];
cubicity31=[0 1/3 2/3 1];

fileID=fopen('gridNoNoise.txt','w');
bigGrid=zeros(600,6);
for i=1:600
   bigGrid(i,:)=[density(randi(9)) angle1(randi(4)) angle2(randi(4)) angle3(randi(4)) cubicity21(randi(4)) cubicity31(randi(4))];
end

for i=1:600
    fprintf(fileID,'%14.10f',bigGrid(i,:));
    fprintf(fileID,'\n');
end