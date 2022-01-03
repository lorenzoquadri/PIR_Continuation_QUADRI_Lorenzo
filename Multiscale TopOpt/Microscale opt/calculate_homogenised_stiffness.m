clear
clc
close all

nelx=9;
nely=9;
nelz=9;
angle1=0;
angle2=0;
angle3=3/8*pi;
transmiLim=0.5;
penal=3;
%% MATERIAL PROPERTIES
E0=1;
Emin=1e-9;
nu=0.3;

%% Evaluation of the initial design
%% PREPARE FINITE ELEMENT ANALYSIS
KE = keGen(nu);
nodenrs = reshape(1:(1+nelx)*(1+nely)*(1+nelz),1+nely,1+nelx,1+nelz);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1, nelx*nely*nelz,1);
edofMat = repmat(edofVec,1,24)+repmat([0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1 3*(nelx+1)*(nely+1)+[0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1]],nelx*nely*nelz,1);
iK = reshape(kron(edofMat,ones(24,1))',576*nelx*nely*nelz,1);
jK = reshape(kron(edofMat,ones(1,24))',576*nelx*nely*nelz,1);

%% PERIODIC BOUNDARY CONDITIONS
e0 = eye(6);
ufixed = zeros(24,6);
U = zeros(3*(nely+1)*(nelx+1)*(nelz+1),6);
alldofs = (1:3*(nely+1)*(nelx+1)*(nelz+1)); 
n1 = [nodenrs(end,[1,end],1),nodenrs(1,[end,1],1),nodenrs(end,[1,end],end),nodenrs(1,[end,1],end)];
d1 = reshape([(3*n1-2);(3*n1-1);3*n1],1,24);
for j = 1:6
    ufixed(4:6,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[nelx;0;0];
    ufixed(7:9,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[nelx;nely;0];
    ufixed(10:12,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[0;nely;0];
    ufixed(13:15,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[0;0;nelz];
    ufixed(16:18,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[nelx;0;nelz];
    ufixed(19:21,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[nelx;nely;nelz];
    ufixed(22:24,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[0;nely;nelz];
end
%% INITIALIZE ITERATION
qe= cell(6,6);
Q=zeros(6,6);
dQ = cell(6,6);

% %Verification using an optimised cell and its symmetric with respect to
% %the plane yz
% x_load=load('micro.mat');
% x=x_load.micro;
% x=fliplr(x); % find symmetric wrt plane yz

% Verification generating a random density distribution, making it central
% symmetric
x = rand(9,9,9);
x=(x+flip(fliplr(flipud(x)),3))/2; %make design central symmetric
% x=fliplr(x); % find symmetric wrt plane yz
xPhys = x;


%obtain side elements used (farthestUsed indicates the first or last
    %element INCLUDED in the determination of the transfer zones)
        
    %top left front corner - top side
    farthestUsedy0x0z0=0;
    xTestedy0x0z0=xPhys(1,farthestUsedy0x0z0+1,1);
    while (xTestedy0x0z0 > transmiLim && farthestUsedy0x0z0<nelx)
        farthestUsedy0x0z0=farthestUsedy0x0z0+1;
        if(farthestUsedy0x0z0<nelx)
            xTestedy0x0z0=xPhys(1,farthestUsedy0x0z0+1,1);
        end
    end
    
    %top left front corner - left side
    farthestUsedx0y0z0=0;
    xTestedx0y0z0=xPhys(farthestUsedx0y0z0+1,1,1);
    while (xTestedx0y0z0 > transmiLim && farthestUsedx0y0z0<nely)
        farthestUsedx0y0z0=farthestUsedx0y0z0+1;
        if(farthestUsedx0y0z0<nely)
            xTestedx0y0z0=xPhys(farthestUsedx0y0z0+1,1,1);
        end
    end
    
    %top left front corner - depth side
    farthestUsedz0x0y0=0;
    xTestedz0x0y0=xPhys(1,1,farthestUsedz0x0y0+1);
    while (xTestedz0x0y0 > transmiLim && farthestUsedz0x0y0<nelz)
        farthestUsedz0x0y0=farthestUsedz0x0y0+1;
        if(farthestUsedz0x0y0<nelz)
            xTestedz0x0y0=xPhys(1,1,farthestUsedz0x0y0+1);
        end
    end
    
    %top right front corner - top side
    farthestUsedy0xMz0=nelx+1;
    xTestedy0xMz0=xPhys(1,farthestUsedy0xMz0-1,1);
    while (xTestedy0xMz0 > transmiLim && farthestUsedy0xMz0>1)
        farthestUsedy0xMz0=farthestUsedy0xMz0-1;
        if(farthestUsedy0xMz0>1)
            xTestedy0xMz0=xPhys(1,farthestUsedy0xMz0-1,1);
        end
    end
    
    %bottom left front corner - left side
    farthestUsedx0yMz0=nely+1;
    xTestedx0yMz0=xPhys(farthestUsedx0yMz0-1,1,1);
    while (xTestedx0yMz0 > transmiLim && farthestUsedx0yMz0>1)
        farthestUsedx0yMz0=farthestUsedx0yMz0-1;
        if(farthestUsedx0yMz0>1)
            xTestedx0yMz0=xPhys(farthestUsedx0yMz0-1,1,1);
        end
    end
    
    %top left rear corner - depth side
    farthestUsedzMx0y0=nelz+1;
    xTestedzMx0y0=xPhys(1,1,farthestUsedzMx0y0-1);
    while (xTestedzMx0y0 > transmiLim && farthestUsedzMx0y0>1)
        farthestUsedzMx0y0=farthestUsedzMx0y0-1;
        if(farthestUsedzMx0y0>1)
            xTestedzMx0y0=xPhys(1,1,farthestUsedzMx0y0-1);
        end
    end
    
    %Make all transmission lengths equal
    lengthtop0=farthestUsedy0x0z0;
    lengthtopM=nelx-farthestUsedy0xMz0+1;
    lengthtop=min([lengthtop0, lengthtopM]);
    farthestUsedy0x0z0=round(lengthtop);
    farthestUsedy0xMz0=nelx-round(lengthtop)+1;
     
    lengthleft0=farthestUsedx0y0z0;
    lengthleftM=nely-farthestUsedx0yMz0+1;
    lengthleft=min([lengthleft0,lengthleftM]);
    farthestUsedx0y0z0=round(lengthleft);
    farthestUsedx0yMz0=nely-round(lengthleft)+1;
    
    lengthdepth0=farthestUsedz0x0y0;
    lengthdepthM=nelz-farthestUsedzMx0y0+1;
    lengthdepth=min([lengthdepth0,lengthdepthM]);
    farthestUsedz0x0y0=round(lengthdepth);
    farthestUsedzMx0y0=nelz-round(lengthdepth)+1;
    
    %make sure no double count
    if (farthestUsedx0y0z0 >= farthestUsedx0yMz0)
        farthestUsedx0y0z0=floor((nely)/2);
        farthestUsedx0yMz0=ceil((nely)/2)+1;
    end
    if (farthestUsedy0x0z0 >= farthestUsedy0xMz0)
        farthestUsedy0x0z0=floor((nelx)/2);
        farthestUsedy0xMz0=ceil((nelx)/2)+1;
    end
    if (farthestUsedz0x0y0 >= farthestUsedzMx0y0)
        farthestUsedz0x0y0=floor((nelz)/2);
        farthestUsedzMx0y0=ceil((nelz)/2)+1;
    end
    
    % calculation of the nodes associated to the elements corresponding to
    % the transmission zones
    
    %determination if even or odd number of elements
    if(mod(nely,2)==0)
        nfarthestUsedx0y0z0=farthestUsedx0y0z0;
        nfarthestUsedx0yMz0=farthestUsedx0yMz0+1;
    else
        nfarthestUsedx0y0z0=farthestUsedx0y0z0+1;
        nfarthestUsedx0yMz0=farthestUsedx0yMz0;
    end
    
    if(mod(nelx,2)==0)
        nfarthestUsedy0x0z0=farthestUsedy0x0z0;
        nfarthestUsedy0xMz0=farthestUsedy0xMz0+1;
    else
        nfarthestUsedy0x0z0=farthestUsedy0x0z0+1;
        nfarthestUsedy0xMz0=farthestUsedy0xMz0;
    end
    
    if(mod(nelz,2)==0)
        nfarthestUsedz0x0y0=farthestUsedz0x0y0;
        nfarthestUsedzMx0y0=farthestUsedzMx0y0+1;
    else
        nfarthestUsedz0x0y0=farthestUsedz0x0y0+1;
        nfarthestUsedzMx0y0=farthestUsedzMx0y0;
    end
    
    %moved n, d and wfixed here, to adapt to size of used side of cell

% %     one matrix next to the other w/o reshaping
%     %n3 = [nodenrs(2:farthestUsedx0y0,1)',nodenrs(farthestUsedx0y1+1:farthestUsedx0y2+1,1)',nodenrs(farthestUsedx0y3+1:farthestUsedx0y4+1,1)',nodenrs(farthestUsedx0yM+1:nely,1)',nodenrs(end,2:farthestUsedy0x0),nodenrs(end,farthestUsedy0x1+1:farthestUsedy0x2+1),nodenrs(end,farthestUsedy0x3+1:farthestUsedy0x4+1),nodenrs(end,farthestUsedy0xM+1:nelx)];
%     n3_matrix = [nodenrs(1:farthestUsedx0y0z0,2:farthestUsedy0x0z0,1),nodenrs(farthestUsedx0yMz0:nely+1,2:farthestUsedy0x0z0,1),nodenrs(farthestUsedx0yMz0:nely+1,farthestUsedy0xMz0:nelx,1),nodenrs(1:farthestUsedx0y0z0,farthestUsedy0xMz0:nelx,1),... %front
%         permute(nodenrs(2:farthestUsedx0y0z0,1,1:farthestUsedz0x0y0),[3,1,2]),permute(nodenrs(farthestUsedx0yMz0:nely,1,1:farthestUsedz0x0y0),[3,1,2]),permute(nodenrs(farthestUsedx0yMz0:nely,1,farthestUsedzMx0y0:nelz+1),[3,1,2]),permute(nodenrs(2:farthestUsedx0y0z0,1,farthestUsedzMx0y0:nelz+1),[3,1,2]),... %left
%         permute(nodenrs(nely+1,1:farthestUsedy0x0z0,2:farthestUsedz0x0y0),[2,3,1]),permute(nodenrs(nely+1,farthestUsedy0xMz0:nelx+1,2:farthestUsedz0x0y0),[2,3,1]),permute(nodenrs(nely+1,farthestUsedy0xMz0:nelx+1,farthestUsedzMx0y0:nelz),[2,3,1]),permute(nodenrs(nely+1,1:farthestUsedy0x0z0,farthestUsedzMx0y0:nelz),[2,3,1])]; %bottom
%     n3=reshape(n3_matrix,size(n3_matrix,1)*size(n3_matrix,2),1);
%     d3 = reshape([(3*n3-2);(3*n3-1);3*n3],1,3*size(n3,1)); 
%     n4_matrix = [nodenrs(1:farthestUsedx0y0z0,2:farthestUsedy0x0z0,nelz+1),nodenrs(farthestUsedx0yMz0:nely+1,2:farthestUsedy0x0z0,nelz+1),nodenrs(farthestUsedx0yMz0:nely+1,farthestUsedy0xMz0:nelx,nelz+1),nodenrs(1:farthestUsedx0y0z0,farthestUsedy0xMz0:nelx,nelz+1),... %rear
%         permute(nodenrs(2:farthestUsedx0y0z0,nelx+1,1:farthestUsedz0x0y0),[3,1,2]),permute(nodenrs(farthestUsedx0yMz0:nely,nelx+1,1:farthestUsedz0x0y0),[3,1,2]),permute(nodenrs(farthestUsedx0yMz0:nely,nelx+1,farthestUsedzMx0y0:nelz+1),[3,1,2]),permute(nodenrs(2:farthestUsedx0y0z0,nelx+1,farthestUsedzMx0y0:nelz+1),[3,1,2]),... % right
%         permute(nodenrs(1,1:farthestUsedy0x0z0,2:farthestUsedz0x0y0),[2,3,1]),permute(nodenrs(1,farthestUsedy0xMz0:nelx+1,2:farthestUsedz0x0y0),[2,3,1]),permute(nodenrs(1,farthestUsedy0xMz0:nelx+1,farthestUsedzMx0y0:nelz),[2,3,1]),permute(nodenrs(1,1:farthestUsedy0x0z0,farthestUsedzMx0y0:nelz),[2,3,1])]; %top
%     n4=reshape(n4_matrix,size(n4_matrix,1)*size(n4_matrix,2),1);
%     d4 = reshape([(3*n4-2);(3*n4-1);3*n4],1,3*size(n4,1)); 
%     d2 = setdiff(alldofs,[d1,d3,d4]);
%     
    
    n3 = [reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,1,1:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,1,1:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,1,nfarthestUsedzMx0y0:nelz+1),[3,1,2]),[],1);reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,1,nfarthestUsedzMx0y0:nelz+1),[3,1,2]),[],1);... %left
        reshape(permute(nodenrs(nely+1,1:nfarthestUsedy0x0z0,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,nfarthestUsedy0xMz0:nelx+1,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,nfarthestUsedy0xMz0:nelx+1,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,1:nfarthestUsedy0x0z0,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);... %bottom
        reshape(nodenrs(1:nfarthestUsedx0y0z0,2:nfarthestUsedy0x0z0,1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely+1,2:nfarthestUsedy0x0z0,1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely+1,nfarthestUsedy0xMz0:nelx,1),[],1);reshape(nodenrs(1:nfarthestUsedx0y0z0,nfarthestUsedy0xMz0:nelx,1),[],1)]'; %front
    d3 = reshape([(3*n3-2);(3*n3-1);3*n3],1,3*size(n3,2)); 
    n4 = [reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,nelx+1,1:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,nelx+1,1:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,nelx+1,nfarthestUsedzMx0y0:nelz+1),[3,1,2]),[],1);reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,nelx+1,nfarthestUsedzMx0y0:nelz+1),[3,1,2]),[],1);... % right
        reshape(permute(nodenrs(1,1:nfarthestUsedy0x0z0,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(1,nfarthestUsedy0xMz0:nelx+1,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(1,nfarthestUsedy0xMz0:nelx+1,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);reshape(permute(nodenrs(1,1:nfarthestUsedy0x0z0,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);... %top
        reshape(nodenrs(1:nfarthestUsedx0y0z0,2:nfarthestUsedy0x0z0,nelz+1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely+1,2:nfarthestUsedy0x0z0,nelz+1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely+1,nfarthestUsedy0xMz0:nelx,nelz+1),[],1);reshape(nodenrs(1:nfarthestUsedx0y0z0,nfarthestUsedy0xMz0:nelx,nelz+1),[],1)]'; %rear
    d4 = reshape([(3*n4-2);(3*n4-1);3*n4],1,3*size(n4,2)); 
    d2 = setdiff(alldofs,[d1,d3,d4]);
       
    % wfixed = [repmat(ufixed(4:6,:),(farthestUsedy0x0z0+nelx+1-farthestUsedy0xMz0)*(farthestUsedz0x0y0+nelz+1-farthestUsedzMx0y0)-1,1); repmat(ufixed(10:12,:),(farthestUsedx0y0z0+nely+1-farthestUsedx0yMz0)*(farthestUsedz0x0y0+nelz+1-farthestUsedzMx0y0)-1,1); repmat(ufixed(13:15,:),(farthestUsedx0y0z0+nely+1-farthestUsedx0yMz0)*(farthestUsedz0x0y0+nelx+1-farthestUsedy0xMz0)-1,1)];
    % wfixed = [repmat(ufixed(4:6,:),(farthestUsedy0x0z0-1+nelx+1-farthestUsedy0xMz0)*(farthestUsedz0x0y0+nely+2-farthestUsedzMx0y0),1); repmat(ufixed(10:12,:),(farthestUsedx0y0z0+nely+2-farthestUsedx0yMz0)*(farthestUsedz0x0y0-1+nelz+1-farthestUsedzMx0y0),1); repmat(ufixed(13:15,:),(farthestUsedx0y0z0+nely+2-farthestUsedx0yMz0)*(farthestUsedz0x0y0-1+nelx+1-farthestUsedy0xMz0),1)];
    
    wfixed = [repmat(ufixed(4:6,:),(nfarthestUsedx0y0z0-1+nely+1-nfarthestUsedx0yMz0)*(nfarthestUsedz0x0y0+nelz+2-nfarthestUsedzMx0y0),1); repmat(ufixed(10:12,:),(nfarthestUsedy0x0z0+nelx+2-nfarthestUsedy0xMz0)*(nfarthestUsedz0x0y0-1+nelz+1-nfarthestUsedzMx0y0),1); repmat(ufixed(13:15,:),(nfarthestUsedx0y0z0+nely+2-nfarthestUsedx0yMz0)*(nfarthestUsedy0x0z0-1+nelx+1-nfarthestUsedy0xMz0),1)];
        
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),576*nelx*nely*nelz,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    Kr = [K(d2,d2), K(d2,d3)+K(d2,d4);K(d3,d2)+K(d4,d2), K(d3,d3)+K(d3,d4)+K(d4,d3)+K(d4,d4)];
    U(d1,:) = ufixed;
    U([d2,d3],:) = Kr\(-[K(d2,d1);K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4); K(d3,d4)+K(d4,d4)]*wfixed);
    
    U(d4,:) = U(d3,:)+wfixed;
    
       
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    for i = 1:6
        for j = 1:6
            U1 = U(:,i); U2 = U(:,j);
            qe{i,j} = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx,nelz)/(nelx*nely*nelz);
            Q(i,j) = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*qe{i,j})));
            dQ{i,j} = penal*(E0-Emin)*xPhys.^(penal-1).*qe{i,j};
        end
    end
    Q
 figure
display_3D(xPhys);

%% Evaluation of the flipped

%% PREPARE FINITE ELEMENT ANALYSIS
KE = keGen(nu);
nodenrs = reshape(1:(1+nelx)*(1+nely)*(1+nelz),1+nely,1+nelx,1+nelz);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1, nelx*nely*nelz,1);
edofMat = repmat(edofVec,1,24)+repmat([0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1 3*(nelx+1)*(nely+1)+[0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1]],nelx*nely*nelz,1);
iK = reshape(kron(edofMat,ones(24,1))',576*nelx*nely*nelz,1);
jK = reshape(kron(edofMat,ones(1,24))',576*nelx*nely*nelz,1);

%% PERIODIC BOUNDARY CONDITIONS
e0 = eye(6);
ufixed = zeros(24,6);
U = zeros(3*(nely+1)*(nelx+1)*(nelz+1),6);
alldofs = (1:3*(nely+1)*(nelx+1)*(nelz+1)); 
n1 = [nodenrs(end,[1,end],1),nodenrs(1,[end,1],1),nodenrs(end,[1,end],end),nodenrs(1,[end,1],end)];
d1 = reshape([(3*n1-2);(3*n1-1);3*n1],1,24);
for j = 1:6
    ufixed(4:6,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[nelx;0;0];
    ufixed(7:9,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[nelx;nely;0];
    ufixed(10:12,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[0;nely;0];
    ufixed(13:15,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[0;0;nelz];
    ufixed(16:18,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[nelx;0;nelz];
    ufixed(19:21,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[nelx;nely;nelz];
    ufixed(22:24,j) = [e0(1,j),e0(6,j)/2,e0(5,j)/2;e0(6,j)/2,e0(2,j),e0(4,j)/2;e0(5,j)/2,e0(4,j)/2,e0(3,j)]*[0;nely;nelz];
end
%% INITIALIZE ITERATION
qe= cell(6,6);
Q=zeros(6,6);
dQ = cell(6,6);
x=fliplr(x); % find symmetric wrt plane yz
xPhys = x;

%obtain side elements used (farthestUsed indicates the first or last
    %element INCLUDED in the determination of the transfer zones)
        
    %top left front corner - top side
    farthestUsedy0x0z0=0;
    xTestedy0x0z0=xPhys(1,farthestUsedy0x0z0+1,1);
    while (xTestedy0x0z0 > transmiLim && farthestUsedy0x0z0<nelx)
        farthestUsedy0x0z0=farthestUsedy0x0z0+1;
        if(farthestUsedy0x0z0<nelx)
            xTestedy0x0z0=xPhys(1,farthestUsedy0x0z0+1,1);
        end
    end
    
    %top left front corner - left side
    farthestUsedx0y0z0=0;
    xTestedx0y0z0=xPhys(farthestUsedx0y0z0+1,1,1);
    while (xTestedx0y0z0 > transmiLim && farthestUsedx0y0z0<nely)
        farthestUsedx0y0z0=farthestUsedx0y0z0+1;
        if(farthestUsedx0y0z0<nely)
            xTestedx0y0z0=xPhys(farthestUsedx0y0z0+1,1,1);
        end
    end
    
    %top left front corner - depth side
    farthestUsedz0x0y0=0;
    xTestedz0x0y0=xPhys(1,1,farthestUsedz0x0y0+1);
    while (xTestedz0x0y0 > transmiLim && farthestUsedz0x0y0<nelz)
        farthestUsedz0x0y0=farthestUsedz0x0y0+1;
        if(farthestUsedz0x0y0<nelz)
            xTestedz0x0y0=xPhys(1,1,farthestUsedz0x0y0+1);
        end
    end
    
    %top right front corner - top side
    farthestUsedy0xMz0=nelx+1;
    xTestedy0xMz0=xPhys(1,farthestUsedy0xMz0-1,1);
    while (xTestedy0xMz0 > transmiLim && farthestUsedy0xMz0>1)
        farthestUsedy0xMz0=farthestUsedy0xMz0-1;
        if(farthestUsedy0xMz0>1)
            xTestedy0xMz0=xPhys(1,farthestUsedy0xMz0-1,1);
        end
    end
    
    %bottom left front corner - left side
    farthestUsedx0yMz0=nely+1;
    xTestedx0yMz0=xPhys(farthestUsedx0yMz0-1,1,1);
    while (xTestedx0yMz0 > transmiLim && farthestUsedx0yMz0>1)
        farthestUsedx0yMz0=farthestUsedx0yMz0-1;
        if(farthestUsedx0yMz0>1)
            xTestedx0yMz0=xPhys(farthestUsedx0yMz0-1,1,1);
        end
    end
    
    %top left rear corner - depth side
    farthestUsedzMx0y0=nelz+1;
    xTestedzMx0y0=xPhys(1,1,farthestUsedzMx0y0-1);
    while (xTestedzMx0y0 > transmiLim && farthestUsedzMx0y0>1)
        farthestUsedzMx0y0=farthestUsedzMx0y0-1;
        if(farthestUsedzMx0y0>1)
            xTestedzMx0y0=xPhys(1,1,farthestUsedzMx0y0-1);
        end
    end
    
    %Make all transmission lengths equal
    lengthtop0=farthestUsedy0x0z0;
    lengthtopM=nelx-farthestUsedy0xMz0+1;
    lengthtop=min([lengthtop0, lengthtopM]);
    farthestUsedy0x0z0=round(lengthtop);
    farthestUsedy0xMz0=nelx-round(lengthtop)+1;
     
    lengthleft0=farthestUsedx0y0z0;
    lengthleftM=nely-farthestUsedx0yMz0+1;
    lengthleft=min([lengthleft0,lengthleftM]);
    farthestUsedx0y0z0=round(lengthleft);
    farthestUsedx0yMz0=nely-round(lengthleft)+1;
    
    lengthdepth0=farthestUsedz0x0y0;
    lengthdepthM=nelz-farthestUsedzMx0y0+1;
    lengthdepth=min([lengthdepth0,lengthdepthM]);
    farthestUsedz0x0y0=round(lengthdepth);
    farthestUsedzMx0y0=nelz-round(lengthdepth)+1;
    
    %make sure no double count
    if (farthestUsedx0y0z0 >= farthestUsedx0yMz0)
        farthestUsedx0y0z0=floor((nely)/2);
        farthestUsedx0yMz0=ceil((nely)/2)+1;
    end
    if (farthestUsedy0x0z0 >= farthestUsedy0xMz0)
        farthestUsedy0x0z0=floor((nelx)/2);
        farthestUsedy0xMz0=ceil((nelx)/2)+1;
    end
    if (farthestUsedz0x0y0 >= farthestUsedzMx0y0)
        farthestUsedz0x0y0=floor((nelz)/2);
        farthestUsedzMx0y0=ceil((nelz)/2)+1;
    end
    
    % calculation of the nodes associated to the elements corresponding to
    % the transmission zones
    
    %determination if even or odd number of elements
    if(mod(nely,2)==0)
        nfarthestUsedx0y0z0=farthestUsedx0y0z0;
        nfarthestUsedx0yMz0=farthestUsedx0yMz0+1;
    else
        nfarthestUsedx0y0z0=farthestUsedx0y0z0+1;
        nfarthestUsedx0yMz0=farthestUsedx0yMz0;
    end
    
    if(mod(nelx,2)==0)
        nfarthestUsedy0x0z0=farthestUsedy0x0z0;
        nfarthestUsedy0xMz0=farthestUsedy0xMz0+1;
    else
        nfarthestUsedy0x0z0=farthestUsedy0x0z0+1;
        nfarthestUsedy0xMz0=farthestUsedy0xMz0;
    end
    
    if(mod(nelz,2)==0)
        nfarthestUsedz0x0y0=farthestUsedz0x0y0;
        nfarthestUsedzMx0y0=farthestUsedzMx0y0+1;
    else
        nfarthestUsedz0x0y0=farthestUsedz0x0y0+1;
        nfarthestUsedzMx0y0=farthestUsedzMx0y0;
    end
    
    %moved n, d and wfixed here, to adapt to size of used side of cell

% %     one matrix next to the other w/o reshaping
%     %n3 = [nodenrs(2:farthestUsedx0y0,1)',nodenrs(farthestUsedx0y1+1:farthestUsedx0y2+1,1)',nodenrs(farthestUsedx0y3+1:farthestUsedx0y4+1,1)',nodenrs(farthestUsedx0yM+1:nely,1)',nodenrs(end,2:farthestUsedy0x0),nodenrs(end,farthestUsedy0x1+1:farthestUsedy0x2+1),nodenrs(end,farthestUsedy0x3+1:farthestUsedy0x4+1),nodenrs(end,farthestUsedy0xM+1:nelx)];
%     n3_matrix = [nodenrs(1:farthestUsedx0y0z0,2:farthestUsedy0x0z0,1),nodenrs(farthestUsedx0yMz0:nely+1,2:farthestUsedy0x0z0,1),nodenrs(farthestUsedx0yMz0:nely+1,farthestUsedy0xMz0:nelx,1),nodenrs(1:farthestUsedx0y0z0,farthestUsedy0xMz0:nelx,1),... %front
%         permute(nodenrs(2:farthestUsedx0y0z0,1,1:farthestUsedz0x0y0),[3,1,2]),permute(nodenrs(farthestUsedx0yMz0:nely,1,1:farthestUsedz0x0y0),[3,1,2]),permute(nodenrs(farthestUsedx0yMz0:nely,1,farthestUsedzMx0y0:nelz+1),[3,1,2]),permute(nodenrs(2:farthestUsedx0y0z0,1,farthestUsedzMx0y0:nelz+1),[3,1,2]),... %left
%         permute(nodenrs(nely+1,1:farthestUsedy0x0z0,2:farthestUsedz0x0y0),[2,3,1]),permute(nodenrs(nely+1,farthestUsedy0xMz0:nelx+1,2:farthestUsedz0x0y0),[2,3,1]),permute(nodenrs(nely+1,farthestUsedy0xMz0:nelx+1,farthestUsedzMx0y0:nelz),[2,3,1]),permute(nodenrs(nely+1,1:farthestUsedy0x0z0,farthestUsedzMx0y0:nelz),[2,3,1])]; %bottom
%     n3=reshape(n3_matrix,size(n3_matrix,1)*size(n3_matrix,2),1);
%     d3 = reshape([(3*n3-2);(3*n3-1);3*n3],1,3*size(n3,1)); 
%     n4_matrix = [nodenrs(1:farthestUsedx0y0z0,2:farthestUsedy0x0z0,nelz+1),nodenrs(farthestUsedx0yMz0:nely+1,2:farthestUsedy0x0z0,nelz+1),nodenrs(farthestUsedx0yMz0:nely+1,farthestUsedy0xMz0:nelx,nelz+1),nodenrs(1:farthestUsedx0y0z0,farthestUsedy0xMz0:nelx,nelz+1),... %rear
%         permute(nodenrs(2:farthestUsedx0y0z0,nelx+1,1:farthestUsedz0x0y0),[3,1,2]),permute(nodenrs(farthestUsedx0yMz0:nely,nelx+1,1:farthestUsedz0x0y0),[3,1,2]),permute(nodenrs(farthestUsedx0yMz0:nely,nelx+1,farthestUsedzMx0y0:nelz+1),[3,1,2]),permute(nodenrs(2:farthestUsedx0y0z0,nelx+1,farthestUsedzMx0y0:nelz+1),[3,1,2]),... % right
%         permute(nodenrs(1,1:farthestUsedy0x0z0,2:farthestUsedz0x0y0),[2,3,1]),permute(nodenrs(1,farthestUsedy0xMz0:nelx+1,2:farthestUsedz0x0y0),[2,3,1]),permute(nodenrs(1,farthestUsedy0xMz0:nelx+1,farthestUsedzMx0y0:nelz),[2,3,1]),permute(nodenrs(1,1:farthestUsedy0x0z0,farthestUsedzMx0y0:nelz),[2,3,1])]; %top
%     n4=reshape(n4_matrix,size(n4_matrix,1)*size(n4_matrix,2),1);
%     d4 = reshape([(3*n4-2);(3*n4-1);3*n4],1,3*size(n4,1)); 
%     d2 = setdiff(alldofs,[d1,d3,d4]);
%     
    
    n3 = [reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,1,1:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,1,1:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,1,nfarthestUsedzMx0y0:nelz+1),[3,1,2]),[],1);reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,1,nfarthestUsedzMx0y0:nelz+1),[3,1,2]),[],1);... %left
        reshape(permute(nodenrs(nely+1,1:nfarthestUsedy0x0z0,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,nfarthestUsedy0xMz0:nelx+1,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,nfarthestUsedy0xMz0:nelx+1,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,1:nfarthestUsedy0x0z0,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);... %bottom
        reshape(nodenrs(1:nfarthestUsedx0y0z0,2:nfarthestUsedy0x0z0,1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely+1,2:nfarthestUsedy0x0z0,1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely+1,nfarthestUsedy0xMz0:nelx,1),[],1);reshape(nodenrs(1:nfarthestUsedx0y0z0,nfarthestUsedy0xMz0:nelx,1),[],1)]'; %front
    d3 = reshape([(3*n3-2);(3*n3-1);3*n3],1,3*size(n3,2)); 
    n4 = [reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,nelx+1,1:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,nelx+1,1:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,nelx+1,nfarthestUsedzMx0y0:nelz+1),[3,1,2]),[],1);reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,nelx+1,nfarthestUsedzMx0y0:nelz+1),[3,1,2]),[],1);... % right
        reshape(permute(nodenrs(1,1:nfarthestUsedy0x0z0,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(1,nfarthestUsedy0xMz0:nelx+1,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(1,nfarthestUsedy0xMz0:nelx+1,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);reshape(permute(nodenrs(1,1:nfarthestUsedy0x0z0,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);... %top
        reshape(nodenrs(1:nfarthestUsedx0y0z0,2:nfarthestUsedy0x0z0,nelz+1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely+1,2:nfarthestUsedy0x0z0,nelz+1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely+1,nfarthestUsedy0xMz0:nelx,nelz+1),[],1);reshape(nodenrs(1:nfarthestUsedx0y0z0,nfarthestUsedy0xMz0:nelx,nelz+1),[],1)]'; %rear
    d4 = reshape([(3*n4-2);(3*n4-1);3*n4],1,3*size(n4,2)); 
    d2 = setdiff(alldofs,[d1,d3,d4]);
       
    % wfixed = [repmat(ufixed(4:6,:),(farthestUsedy0x0z0+nelx+1-farthestUsedy0xMz0)*(farthestUsedz0x0y0+nelz+1-farthestUsedzMx0y0)-1,1); repmat(ufixed(10:12,:),(farthestUsedx0y0z0+nely+1-farthestUsedx0yMz0)*(farthestUsedz0x0y0+nelz+1-farthestUsedzMx0y0)-1,1); repmat(ufixed(13:15,:),(farthestUsedx0y0z0+nely+1-farthestUsedx0yMz0)*(farthestUsedz0x0y0+nelx+1-farthestUsedy0xMz0)-1,1)];
    % wfixed = [repmat(ufixed(4:6,:),(farthestUsedy0x0z0-1+nelx+1-farthestUsedy0xMz0)*(farthestUsedz0x0y0+nely+2-farthestUsedzMx0y0),1); repmat(ufixed(10:12,:),(farthestUsedx0y0z0+nely+2-farthestUsedx0yMz0)*(farthestUsedz0x0y0-1+nelz+1-farthestUsedzMx0y0),1); repmat(ufixed(13:15,:),(farthestUsedx0y0z0+nely+2-farthestUsedx0yMz0)*(farthestUsedz0x0y0-1+nelx+1-farthestUsedy0xMz0),1)];
    
    wfixed = [repmat(ufixed(4:6,:),(nfarthestUsedx0y0z0-1+nely+1-nfarthestUsedx0yMz0)*(nfarthestUsedz0x0y0+nelz+2-nfarthestUsedzMx0y0),1); repmat(ufixed(10:12,:),(nfarthestUsedy0x0z0+nelx+2-nfarthestUsedy0xMz0)*(nfarthestUsedz0x0y0-1+nelz+1-nfarthestUsedzMx0y0),1); repmat(ufixed(13:15,:),(nfarthestUsedx0y0z0+nely+2-nfarthestUsedx0yMz0)*(nfarthestUsedy0x0z0-1+nelx+1-nfarthestUsedy0xMz0),1)];
        
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),576*nelx*nely*nelz,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    Kr = [K(d2,d2), K(d2,d3)+K(d2,d4);K(d3,d2)+K(d4,d2), K(d3,d3)+K(d3,d4)+K(d4,d3)+K(d4,d4)];
    U(d1,:) = ufixed;
    U([d2,d3],:) = Kr\(-[K(d2,d1);K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4); K(d3,d4)+K(d4,d4)]*wfixed);
    
    U(d4,:) = U(d3,:)+wfixed;
    
       
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    for i = 1:6
        for j = 1:6
            U1 = U(:,i); U2 = U(:,j);
            qe{i,j} = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx,nelz)/(nelx*nely*nelz);
            Q(i,j) = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*qe{i,j})));
            dQ{i,j} = penal*(E0-Emin)*xPhys.^(penal-1).*qe{i,j};
        end
    end
    Q
    figure
display_3D(xPhys);