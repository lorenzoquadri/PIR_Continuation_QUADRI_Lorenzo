clear
clc 
close all
K_optimised=load('K_optimised.mat');
K=load('K.mat');
loop=1;

for i=1:3000
    for j=1:3000
        if abs(K_optimised.K(i,j)-K.K(i,j))>1e-5 && loop==1
            display('error');
            i
            j
            loop=2;
        end
    end
end