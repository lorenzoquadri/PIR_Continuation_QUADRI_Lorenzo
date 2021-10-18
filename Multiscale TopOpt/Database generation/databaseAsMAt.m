clear;clc;
discr=4;

dbfilename='FinalDB/databaseGridNoNoisetot.txt';
dbfileID=fopen(dbfilename);
dbdata = fscanf(dbfileID,'%f',[27 inf]);

dbMat = reshape(dbdata,27,discr,discr,discr,discr,discr,discr);
angle1Length=size(dbMat,6);
angles1=((1:angle1Length)-1)/(angle1Length-1);
angle2Length=size(dbMat,5);
angles2=((1:angle2Length)-1)/(angle2Length-1);
angle3Length=size(dbMat,4);
angles3=((1:angle3Length)-1)/(angle3Length-1);
cub21Length=size(dbMat,3);
cubs21=((1:cub21Length)-1)/(cub21Length-1);
cub31Length=size(dbMat,2);
cubs31=((1:cub31Length)-1)/(cub31Length-1);


%correct cases dens=0 or 1
%lastPoint= [1.09846159394850,0.329428548373026,5.85895174775321e-08,1.09846159394851,5.85895196977880e-08,0.384544699731802];
lastPoint=[1.26269304782024,0.511903774455586,0.511903774455586,-7.88355407931949e-18,-4.40611316683323e-15,-4.45474889119004e-15, 1.26269304782024,0.511903774455586,-4.44138783005208e-15,-1.06460934760249e-18,-4.49859801837204e-15, 1.26269304782024,-4.45523974552029e-15,-4.33213150754008e-15,2.51205963786631e-17, 0.374171338573729,-2.22992623662511e-15,-2.24381541585689e-15, 0.374171338573729,-2.19271665299316e-15, 0.374171338573729];

dbMatFull=zeros(27,discr,discr,discr,discr,discr,discr+2);
dbMatFull(:,:,:,:,:,:,2:end-1)=dbMat;
% for k=1:angle3Length
%     for j=1:angle2Length
%         for i=1:angle1Length
%             for l=1:cub21Length
%                 dbMatFull(:,l,:,k,j,i,1) = [zeros(1,cub21Length);angles1(i)*ones(1,cub21Length);angles2(j)*ones(1,cub21Length);angles3(k)*ones(1,cub21Length);cubs21(:)';zeros(1,cub21Length);zeros(21,cub21Length)];
%                 dbMatFull(:,l,:,k,j,i,end) = [ones(1,cub21Length);angles1(i)*ones(1,cub21Length);angles2(j)*ones(1,cub21Length);angles3(k)*ones(1,cub21Length);cubs21(:)';zeros(1,cub21Length);repmat(lastPoint(:),1,cub21Length)];
%             end
%             for l=1:cub31Length
%                 dbMatFull(:,:,l,k,j,i,1) = dbMatFull(:,:,l,k,j,i,1)+[zeros(1,cub31Length);angles1(i)*ones(1,cub31Length);angles2(j)*ones(1,cub31Length);angles3(k)*ones(1,cub31Length);zeros(1,cub31Length);cubs31(:)';zeros(21,cub31Length)];
%                 dbMatFull(:,:,l,k,j,i,end) = dbMatFull(:,:,l,k,j,i,end)+[ones(1,cub31Length);angles1(i)*ones(1,cub31Length);angles2(j)*ones(1,cub31Length);angles3(k)*ones(1,cub31Length);zeros(1,cub31Length);cubs31(:)';repmat(lastPoint(:),1,cub31Length)];
%             end
%             
%         end
%     end
% end

for i=1:angle1Length
    for j=1:angle2Length
        for k=1:angle3Length
            for l=1:cub21Length
                for m=1:cub31Length
                    
                    dbMatFull(:,m,l,k,j,i,1) = [0;angles1(i);angles2(j);angles3(k);cubs21(l);cubs31(m);zeros(21,1)];
                    dbMatFull(:,m,l,k,j,i,end) = [1;angles1(i);angles2(j);angles3(k);cubs21(l);cubs31(m);lastPoint(:)];
                    
                end
            end
            
        end
    end
end

M=dbMatFull;

%% Expansion of database through symmetry of 45 degrees

%put angle in [0,pi/2] instead of [0,pi/4]
M(2,:,:,:,:,:,:,:)=M(2,:,:,:,:,:,:,:)/2;
M(3,:,:,:,:,:,:,:)=M(3,:,:,:,:,:,:,:)/2;
M(4,:,:,:,:,:,:,:)=M(4,:,:,:,:,:,:,:)/2;

%add symmetric of each point wrt 45 degrees
angle1Length=size(M,6);


Msym1=zeros(size(M,1),size(M,2),size(M,3),size(M,4),size(M,5),size(M,6)-1,size(M,7));


for i=1:size(M,7) %density
    for j=2:size(M,6) %angle1
        for k=1:size(M,5) %angle2
            for l=1:size(M,4) %angle3
                for m=1:size(M,3) %cub21
                    for n=1:size(M,2) %cub31
                        % flip wrt 45 of angle1 around x
                        origPoint1=M(:,n,m,l,k,angle1Length-j+1,i);
                        symPoint1=origPoint1;
                        symPoint1(2)=1-origPoint1(2); %substiution of angle1
                        
                        %inversions of the terms of the elastic tensor
                        symPoint1(13)=origPoint1(18);
                        symPoint1(18)=origPoint1(13);
                        symPoint1(8)=origPoint1(9);
                        symPoint1(9)=origPoint1(8);
                        symPoint1(15)=origPoint1(19);
                        symPoint1(16)=origPoint1(20);
                        symPoint1(17)=origPoint1(21);
                        symPoint1(19)=origPoint1(15);
                        symPoint1(20)=origPoint1(16);
                        symPoint1(21)=origPoint1(17);
                        
                        Msym1(:,n,m,l,k,j-1,i)=symPoint1(:);
                    end
                end
            end   
        end
    end
end
M=cat(6,M,Msym1);

angle2Length=size(M,5);
Msym2=zeros(size(M,1),size(M,2),size(M,3),size(M,4),size(M,5)-1,size(M,6),size(M,7));


for i=1:size(M,7) %density
    for j=1:size(M,6) %angle1
        for k=2:size(M,5) %angle2
            for l=1:size(M,4) %angle3
                for m=1:size(M,3) %cub21
                    for n=1:size(M,2) %cub31
                        % flip wrt 45 of angle2 around y
                        origPoint2=M(:,n,m,l,angle2Length-k+1,j,i);
                        symPoint2=origPoint2;
                        symPoint2(3)=1-origPoint2(3); %substiution of angle2
                        
                        %inversions of the terms of the elastic tensor
                        symPoint2(7)=origPoint2(18);
                        symPoint2(18)=origPoint2(7);
                        symPoint2(8)=origPoint2(14);
                        symPoint2(9)=origPoint2(19);
                        symPoint2(11)=origPoint2(20);
                        symPoint2(12)=origPoint2(21);
                        symPoint2(14)=origPoint2(8);
                        symPoint2(19)=origPoint2(10);
                        symPoint2(20)=origPoint2(11);
                        symPoint2(21)=origPoint2(12);
                        
                        Msym2(:,n,m,l,k-1,j,i)=symPoint2(:);
                    end
                end
            end   
        end
    end
end
M=cat(5,M,Msym2);

angle3Length=size(M,4);
Msym3=zeros(size(M,1),size(M,2),size(M,3),size(M,4)-1,size(M,5),size(M,6),size(M,7));
for i=1:size(M,7) %density
    for j=1:size(M,6) %angle1
        for k=1:size(M,5) %angle2
            for l=2:size(M,4) %angle3
                for m=1:size(M,3) %cub21
                    for n=1:size(M,2) %cub31
                        
                        % flip wrt 45 of angle3 around z
                        origPoint3=M(:,n,m,angle3Length-l+1,k,j,i);
                        symPoint3=origPoint3;
                        symPoint3(4)=1-origPoint3(4); %substiution of angle2
                        
                        %inversions of the terms of the elastic tensor
                        symPoint3(7)=origPoint3(13);
                        symPoint3(13)=origPoint3(7);
                        symPoint3(9)=origPoint3(14);
                        symPoint3(9)=origPoint3(15);
                        symPoint3(11)=origPoint3(16);
                        symPoint3(12)=origPoint3(17);
                        symPoint3(14)=origPoint3(9);
                        symPoint3(15)=origPoint3(10);
                        symPoint3(16)=origPoint3(11);
                        symPoint3(17)=origPoint3(12);
                        
                        Msym3(:,n,m,l-1,k,j,i)=symPoint3(:);
                    end
                end
            end   
        end
    end
end

%M=[M;Msym];
M=cat(4,M,Msym3);


%% Expansion of database through rotation of 90 degrees

%put angle in [0,pi] instead of [0,pi/2]
M(2,:,:,:,:,:,:)=M(2,:,:,:,:,:,:)/2;
M(3,:,:,:,:,:,:,:)=M(3,:,:,:,:,:,:,:)/2;
M(4,:,:,:,:,:,:,:)=M(4,:,:,:,:,:,:,:)/2;

%add 90° rotation of each point
angle1Length=size(M,6);
Mrot1=zeros(size(M,1),size(M,2),size(M,3),size(M,4),size(M,5),size(M,6)-1,size(M,7));


for i=1:size(M,7) %density
    for j=2:size(M,6) %angle1
        for k=1:size(M,5) %angle2
            for l=1:size(M,4) %angle3
                for m=1:size(M,3) %cub21
                    for n=1:size(M,2) %cub31
                        % rotation of angle1 around x
                        origPoint1=M(:,n,m,l,k,j,i);
                        Qin_array1=origPoint1(7:end);
                        Qin_tens1=array2matrix(Qin_array1);
                        Qrot_tens1=rotateTensorMatrix_3D(Qin_tens1,pi/2,0,0);
                        Qrot_array1=matrix2array(Qrot_tens1);
                        rotPoint1=[origPoint1(1) 0.5+origPoint1(2) origPoint1(3) origPoint1(4) origPoint1(5) origPoint1(6) Qrot_array1];
                        Mrot1(:,n,m,l,k,j-1,i)=rotPoint1(:);   
                    end
                end
            end
        end
    end
end
M=cat(6,M,Mrot1);

angle2Length=size(M,5);
Mrot2=zeros(size(M,1),size(M,2),size(M,3),size(M,4),size(M,5)-1,size(M,6),size(M,7));
for i=1:size(M,7) %density
    for j=1:size(M,6) %angle1
        for k=2:size(M,5) %angle2
            for l=1:size(M,4) %angle3
                for m=1:size(M,3) %cub21
                    for n=1:size(M,2) %cub31
                        % rotation of angle2 around y
                        origPoint2=M(:,n,m,l,k,j,i);
                        Qin_array2=origPoint2(7:end);
                        Qin_tens2=array2matrix(Qin_array2);
                        Qrot_tens2=rotateTensorMatrix_3D(Qin_tens2,0,pi/2,0);
                        Qrot_array2=matrix2array(Qrot_tens2);
                        rotPoint2=[origPoint2(1) origPoint2(2) 0.5+origPoint2(3) origPoint2(4) origPoint2(5) origPoint2(6) Qrot_array2];
                        Mrot2(:,n,m,l,k-1,j,i)=rotPoint2(:);
                    end
                end
            end
        end
    end
end
M=cat(5,M,Mrot2);

angle3Length=size(M,4);
Mrot3=zeros(size(M,1),size(M,2),size(M,3),size(M,4)-1,size(M,5),size(M,6),size(M,7));
for i=1:size(M,7) %density
    for j=1:size(M,6) %angle1
        for k=1:size(M,5) %angle2
            for l=2:size(M,4) %angle3
                for m=1:size(M,3) %cub21
                    for n=1:size(M,2) %cub31
                        % rotation of angle3 around z
                        origPoint3=M(:,n,m,l,k,j,i);
                        Qin_array3=origPoint3(7:end);
                        Qin_tens3=array2matrix(Qin_array3);
                        Qrot_tens3=rotateTensorMatrix_3D(Qin_tens3,0,0,pi/2);
                        Qrot_array3=matrix2array(Qrot_tens3);
                        rotPoint3=[origPoint3(1) origPoint3(2) origPoint3(3) 0.5+origPoint3(4) origPoint3(5) origPoint3(6) Qrot_array3];
                        Mrot3(:,n,m,l-1,k,j,i)=rotPoint3(:);
                    end
                end
            end
        end
    end
end

M=cat(4,M,Mrot3);


%% Expansion of database through addition of cubicities
%put principal direction in 1 - 2 instead of 1 - mixed
M(5,:,:,:,:,:,:)=M(5,:,:,:,:,:,:)/2; 
M(6,:,:,:,:,:,:)=M(6,:,:,:,:,:,:)/2; 

%add 90° rotation of each point
cub21Length=size(M,3);
Mrotcub21=zeros(size(M,1),size(M,2),size(M,3)-1,size(M,4),size(M,5),size(M,6),size(M,7));

for i=1:size(M,7) %density
    for j=1:size(M,6) %angle1
        for k=1:size(M,5) %angle2
            for l=1:size(M,4) %angle3
                for m=2:size(M,3) %cub21
                    for n=1:size(M,2) %cub31
                                             
                        
                        % rotation of cub21
                        origPoint21=M(:,n,cub21Length-m+1,l,k,j,i);
                        Qin_array21=origPoint21(7:end);
                        Qin_tens21=array2matrix(Qin_array21);
                        Qrot_tens21=rotateTensorMatrix_3D(Qin_tens21,0,0,pi/2);
                        Qrot_array21=matrix2array(Qrot_tens21);
                        rotPoint21=[origPoint21(1) origPoint21(2) origPoint21(3) origPoint21(4) 1-origPoint21(5) origPoint21(6) Qrot_array21];
                        Mrotcub21(:,n,m-1,l,k,j,i)=rotPoint21(:);
                    end
                end
            end
        end
    end
end
M=cat(3,M,Mrotcub21);

cub31Length=size(M,2);
Mrotcub31=zeros(size(M,1),size(M,2)-1,size(M,3),size(M,4),size(M,5),size(M,6),size(M,7));

for i=1:size(M,7) %density
    for j=1:size(M,6) %angle1
        for k=1:size(M,5) %angle2
            for l=1:size(M,4) %angle3
                for m=1:size(M,3) %cub21
                    for n=2:size(M,2) %cub31
                                               
                        % rotation of cub31
                        origPoint31=M(:,cub31Length-n+1,m,l,k,j,i);
                        Qin_array31=origPoint31(7:end);
                        Qin_tens31=array2matrix(Qin_array31);
                        Qrot_tens31=rotateTensorMatrix_3D(Qin_tens31,0,pi/2,0);
                        Qrot_array31=matrix2array(Qrot_tens31);
                        rotPoint31=[origPoint31(1) origPoint31(2) origPoint31(3) origPoint31(4) origPoint31(5) 1-origPoint31(6) Qrot_array31];
                        Mrotcub31(:,n-1,m,l,k,j,i)=rotPoint31(:);
                    end
                end
            end
        end
    end
end


M=cat(2,M,Mrotcub31);

dbMat=M;

save('database4-4-4-4-4-4.mat', 'dbMat');