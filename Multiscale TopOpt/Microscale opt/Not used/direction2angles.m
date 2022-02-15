function [theta1,theta2,theta3]=direction2angles(X,Y)
theta1=asin(Y(3)/sqrt(1-X(3)^2));
theta2=asin(-X(3));
theta3=asin(X(2)/sqrt(1-X(3)^2));