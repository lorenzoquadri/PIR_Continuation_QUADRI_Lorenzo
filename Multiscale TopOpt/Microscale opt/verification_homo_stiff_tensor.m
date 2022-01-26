clear
clc
close all

l=[1:40]; %i
l_prime=[1:40]; %j
rho=zeros(length(l),length(l_prime));

for i=1:length(l)
    for j=1:i
        if (mod(i-j,2)==0)
            rho_now=1-(l_prime(j)^3+6*l_prime(j)^2*(l(i)-l_prime(j))/2)/l(i)^3;
            rho(i,j)=rho_now;
        end            
    end
end

rho=rho+diag(diag(ones(length(rho))));

delta_rho=rho*100-floor(rho*100);

for i=1:length(l)
    for j=1:length(l_prime)
      if(delta_rho(i,j)==0)
          delta_rho(i,j)=1;
      end
    end
end

[row, col] = find(delta_rho == min(delta_rho(:)));

