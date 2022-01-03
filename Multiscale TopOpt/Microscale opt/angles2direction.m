function [X,Y]=angles2direction(theta1,theta2,theta3)
X(3)=-sin(theta2);
X(2)=sin(theta3)*sqrt(1-X(3)^2);
X(1)=sqrt(1-X(2)^2-X(3)^2);
Y(3)=sin(theta1)*sqrt(1-X(3)^2);