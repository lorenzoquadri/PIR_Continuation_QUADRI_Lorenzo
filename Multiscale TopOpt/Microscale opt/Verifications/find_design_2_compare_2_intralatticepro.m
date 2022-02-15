clear
clc
close all

for nel=20:40
    nelx=nel;
    nely=nel;
    nelz=nel;
    for transmLim=1:10
        [~,rho]=initDesMore8tz_3D(nelx,nely,nelz,0.3,10,transmLim);
        rho100=100*rho;
        rhoint=fix(rho100);
        if abs(rho100-rhoint)<0.02 &rho<0.5
            nel
            transmLim
        end
    end
    
end