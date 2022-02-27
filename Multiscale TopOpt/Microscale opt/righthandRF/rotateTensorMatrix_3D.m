%rotate a 6x6 matrix representing a 4th order 3*3*3*3 elastic tensor
function B=rotateTensorMatrix_3D(A,theta1,theta2,theta3)
B=zeros(6,6);
s1=sin(theta1);
c1=cos(theta2);
s2=sin(theta2);
c2=cos(theta2);
s3=sin(theta3);
c3=cos(theta3);
% Mrot=[[c1*c2*c3-s1*s3 -c1*c2*s3-s1*c3 c1*s2];[s1*c2*c3+c1*s3 -s1*c2*s3+c1*c3 s1*s2];[-s2*c3 s2*s3 c2]];

% Mrot=[[c3*c2 c3*s2*s1-s3*c1 c3*s2*c1+s3*s1];[s3*c2 s3*s2*s1+c3*c2 s3*s2*c1-c3*s1];[-s2 c2*s1 c2*c1]];
% Extrinsic rotation whose improper euler angles are theta1,theta2,theta3
% about axes x,y,z and it's the result of combining three rotation
% matrices: R=Rz(theta3),Ry(theta2)Rx(theta1)

Mrot=[[c2*c3 -c2*s3 s2];[c1*s3+c3*s1*s2 c1*c3-s1*s2*s3 -c2*s1];[s1*s3-c1*c3*s2 c3*s1+c1*s2*s3 c1*c2]];
% Intrinsic rotation whose euler angles are theta1,theta2,theta3
% about intrinsic axes z-y'-x'', which corresponds to a chained rotation 
% about extrinsic axis x-y-z. It's the result of combining three rotation
% matrices: R=Rz(theta3),Ry(theta2)Rx(theta1)

tensA=zeros(3,3,3,3);
tensA(1,1,:,:)=[[A(1,1) A(1,6) A(1,5)];[A(1,6) A(1,2) A(1,4)];[A(1,5) A(1,4) A(1,3)]];
tensA(1,2,:,:)=[[A(6,1) A(6,6) A(6,5)];[A(6,6) A(6,2) A(6,4)];[A(6,5) A(6,4) A(6,3)]];
tensA(1,3,:,:)=[[A(5,1) A(5,6) A(5,5)];[A(5,6) A(5,2) A(5,4)];[A(5,5) A(5,4) A(5,3)]];
tensA(2,1,:,:)=[[A(6,1) A(6,6) A(6,5)];[A(6,6) A(6,2) A(6,4)];[A(6,5) A(6,4) A(6,3)]];
tensA(2,2,:,:)=[[A(2,1) A(2,6) A(2,5)];[A(2,6) A(2,2) A(2,4)];[A(2,5) A(2,4) A(2,3)]];
tensA(2,3,:,:)=[[A(4,1) A(4,6) A(4,5)];[A(4,6) A(4,2) A(4,4)];[A(4,5) A(4,4) A(4,3)]];
tensA(3,1,:,:)=[[A(5,1) A(5,6) A(5,5)];[A(5,6) A(5,2) A(5,4)];[A(5,5) A(5,4) A(5,3)]];
tensA(3,2,:,:)=[[A(4,1) A(4,6) A(4,5)];[A(4,6) A(4,2) A(4,4)];[A(4,5) A(4,4) A(4,3)]];
tensA(3,3,:,:)=[[A(3,1) A(3,6) A(3,5)];[A(3,6) A(3,2) A(3,4)];[A(3,5) A(3,4) A(3,3)]];

tensB=zeros(3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for h=1:3
                for m=1:3
                    for n=1:3
                        for p=1:3
                            for q=1:3
                                tensB(i,j,k,h)=tensB(i,j,k,h)+tensA(m,n,p,q)*Mrot(i,m)*Mrot(j,n)*Mrot(k,p)*Mrot(h,q);
                            end
                        end
                    end
                end
            end
        end
    end
end

B(1,1)=tensB(1,1,1,1);
B(1,2)=tensB(1,1,2,2);
B(1,3)=tensB(1,1,3,3);
B(1,4)=tensB(1,1,2,3);
B(1,5)=tensB(1,1,1,3);
B(1,6)=tensB(1,1,1,2);
B(2,1)=tensB(2,2,1,1);
B(2,2)=tensB(2,2,2,2);
B(2,3)=tensB(2,2,3,3);
B(2,4)=tensB(2,2,2,3);
B(2,5)=tensB(2,2,1,3);
B(2,6)=tensB(2,2,1,2);
B(3,1)=tensB(3,3,1,1);
B(3,2)=tensB(3,3,2,2);
B(3,3)=tensB(3,3,3,3);
B(3,4)=tensB(3,3,2,3);
B(3,5)=tensB(3,3,1,3);
B(3,6)=tensB(3,3,1,2);
B(4,1)=tensB(2,3,1,1);
B(4,2)=tensB(2,3,2,2);
B(4,3)=tensB(2,3,3,3);
B(4,4)=tensB(2,3,2,3);
B(4,5)=tensB(2,3,1,3);
B(4,6)=tensB(2,3,1,2);
B(5,1)=tensB(1,3,1,1);
B(5,2)=tensB(1,3,2,2);
B(5,3)=tensB(1,3,3,3);
B(5,4)=tensB(1,3,2,3);
B(5,5)=tensB(1,3,1,3);
B(5,6)=tensB(1,3,1,2);
B(6,1)=tensB(1,2,1,1);
B(6,2)=tensB(1,2,2,2);
B(6,3)=tensB(1,2,3,3);
B(6,4)=tensB(1,2,2,3);
B(6,5)=tensB(1,2,1,3);
B(6,6)=tensB(1,2,1,2);
end



