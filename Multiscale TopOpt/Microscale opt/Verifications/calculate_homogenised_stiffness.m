clear
clc
close all


nel=20;
nelx=nel;
nely=nel;
nelz=nel;
angle1=0;
angle2=0;
angle3=0;
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

[x,rho]=initDesMore8tz_3D(nelx,nely,nelz,0.3,10,4);

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

n3 = [reshape(nodenrs(nely+1,2:nfarthestUsedy0x0z0,1),[],1);reshape(nodenrs(nely+1,nfarthestUsedy0xMz0:nelx,1),[],1);reshape(nodenrs(2:nfarthestUsedx0y0z0,1,1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely,1,1),[],1);reshape(permute(nodenrs(nely+1,1,2:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nely+1,1,nfarthestUsedzMx0y0:nelz),[3,1,2]),[],1)]';
    d3 = reshape([(3*n3-2);(3*n3-1);3*n3],1,3*size(n3,2)); 
    n4 = [reshape(nodenrs(1,2:nfarthestUsedy0x0z0,1),[],1);reshape(nodenrs(1,nfarthestUsedy0xMz0:nelx,1),[],1);reshape(nodenrs(2:nfarthestUsedx0y0z0,nelx+1,1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely,nelx+1,1),[],1);reshape(permute(nodenrs(nely+1,nelx+1,2:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nely+1,nelx+1,nfarthestUsedzMx0y0:nelz),[3,1,2]),[],1)]';
    d4 = reshape([(3*n4-2);(3*n4-1);3*n4],1,3*size(n4,2)); 
    n5 = [reshape(nodenrs(1,2:nfarthestUsedy0x0z0,nelz+1),[],1);reshape(nodenrs(1,nfarthestUsedy0xMz0:nelx,nelz+1),[],1);reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,nelx+1,nelz+1),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,nelx+1,nelz+1),[3,1,2]),[],1);reshape(permute(nodenrs(1,nelx+1,2:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(1,nelx+1,nfarthestUsedzMx0y0:nelz),[3,1,2]),[],1)]';
    d5 = reshape([(3*n5-2);(3*n5-1);3*n5],1,3*size(n5,2)); 
    n6 = [reshape(nodenrs(nely+1,2:nfarthestUsedy0x0z0,nelz+1),[],1);reshape(nodenrs(nely+1,nfarthestUsedy0xMz0:nelx,nelz+1),[],1);reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,1,nelz+1),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,1,nelz+1),[3,1,2]),[],1);reshape(permute(nodenrs(1,1,2:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(1,1,nfarthestUsedzMx0y0:nelz),[3,1,2]),[],1)]';
    d6 = reshape([(3*n6-2);(3*n6-1);3*n6],1,3*size(n6,2)); 
    n7 = [reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,1,2:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,1,2:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,1,nfarthestUsedzMx0y0:nelz),[3,1,2]),[],1);reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,1,nfarthestUsedzMx0y0:nelz),[3,1,2]),[],1);... %left
         reshape(permute(nodenrs(nely+1,2:nfarthestUsedy0x0z0,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,nfarthestUsedy0xMz0:nelx,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,nfarthestUsedy0xMz0:nelx,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);reshape(permute(nodenrs(nely+1,2:nfarthestUsedy0x0z0,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);... %bottom
         reshape(nodenrs(2:nfarthestUsedx0y0z0,2:nfarthestUsedy0x0z0,1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely,2:nfarthestUsedy0x0z0,1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely,nfarthestUsedy0xMz0:nelx,1),[],1);reshape(nodenrs(2:nfarthestUsedx0y0z0,nfarthestUsedy0xMz0:nelx,1),[],1)]'; %front
    d7 = reshape([(3*n7-2);(3*n7-1);3*n7],1,3*size(n7,2)); 
    n8 = [reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,nelx+1,2:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,nelx+1,2:nfarthestUsedz0x0y0),[3,1,2]),[],1);reshape(permute(nodenrs(nfarthestUsedx0yMz0:nely,nelx+1,nfarthestUsedzMx0y0:nelz),[3,1,2]),[],1);reshape(permute(nodenrs(2:nfarthestUsedx0y0z0,nelx+1,nfarthestUsedzMx0y0:nelz),[3,1,2]),[],1);... % right
        reshape(permute(nodenrs(1,2:nfarthestUsedy0x0z0,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(1,nfarthestUsedy0xMz0:nelx,2:nfarthestUsedz0x0y0),[2,3,1]),[],1);reshape(permute(nodenrs(1,nfarthestUsedy0xMz0:nelx,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);reshape(permute(nodenrs(1,2:nfarthestUsedy0x0z0,nfarthestUsedzMx0y0:nelz),[2,3,1]),[],1);... %top
        reshape(nodenrs(2:nfarthestUsedx0y0z0,2:nfarthestUsedy0x0z0,nelz+1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely,2:nfarthestUsedy0x0z0,nelz+1),[],1);reshape(nodenrs(nfarthestUsedx0yMz0:nely,nfarthestUsedy0xMz0:nelx,nelz+1),[],1);reshape(nodenrs(2:nfarthestUsedx0y0z0,nfarthestUsedy0xMz0:nelx,nelz+1),[],1)]'; %rear
    d8 = reshape([(3*n8-2);(3*n8-1);3*n8],1,3*size(n8,2));  
    d2 = setdiff(alldofs,[d1,d3,d4,d5,d6,d7,d8]);

    w1 = [repmat(ufixed(10:12,:),(nfarthestUsedy0x0z0-1)+(nelx+1-nfarthestUsedy0xMz0),1);repmat(ufixed(4:6,:),(nfarthestUsedx0y0z0-1)+(nely+1-nfarthestUsedx0yMz0),1);  repmat(ufixed(4:6,:),(nfarthestUsedz0x0y0-1)+(nelz+1-nfarthestUsedzMx0y0),1)];
    w2 = [repmat(ufixed(22:24,:),(nfarthestUsedy0x0z0-1)+(nelx+1-nfarthestUsedy0xMz0),1);repmat(ufixed(16:18,:),(nfarthestUsedx0y0z0-1)+(nely+1-nfarthestUsedx0yMz0),1);  repmat(ufixed(7:9,:),(nfarthestUsedz0x0y0-1)+(nelz+1-nfarthestUsedzMx0y0),1)];
    w3 = [repmat(ufixed(13:15,:),(nfarthestUsedy0x0z0-1)+(nelx+1-nfarthestUsedy0xMz0),1);repmat(ufixed(13:15,:),(nfarthestUsedx0y0z0-1)+(nely+1-nfarthestUsedx0yMz0),1);  repmat(ufixed(10:12,:),(nfarthestUsedz0x0y0-1)+(nelz+1-nfarthestUsedzMx0y0),1)];
    w4 = [repmat(ufixed(4:6,:),(nfarthestUsedx0y0z0-1+nely+1-nfarthestUsedx0yMz0)*(nfarthestUsedz0x0y0-1+nelz+1-nfarthestUsedzMx0y0),1); repmat(ufixed(10:12,:),(nfarthestUsedy0x0z0-1+nelx+1-nfarthestUsedy0xMz0)*(nfarthestUsedz0x0y0-1+nelz+1-nfarthestUsedzMx0y0),1); repmat(ufixed(13:15,:),(nfarthestUsedx0y0z0-1+nely+1-nfarthestUsedx0yMz0)*(nfarthestUsedy0x0z0-1+nelx+1-nfarthestUsedy0xMz0),1)];
    
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),576*nelx*nely*nelz,1);

    K = sparse(iK,jK,sK); K = (K+K')/2;

    Kr = [K(d2,d2), K(d2,d3)+K(d2,d4)+K(d2,d5)+K(d2,d6), K(d2,d7)+K(d2,d8);
          K(d3,d2)+K(d4,d2)+K(d5,d2)+K(d6,d2), K(d3,d3)+K(d3,d4)+K(d3,d5)+K(d3,d6)+K(d4,d3)+K(d4,d4)+K(d4,d5)+K(d4,d6)+K(d5,d3)+K(d5,d4)+K(d5,d5)+K(d5,d6)+K(d6,d3)+K(d6,d4)+K(d6,d5)+K(d6,d6), K(d3,d7)+K(d3,d8)+K(d4,d7)+K(d4,d8)+K(d5,d7)+K(d5,d8)+K(d6,d7)+K(d6,d8);
          K(d7,d2)+K(d8,d2), K(d7,d3)+K(d8,d3)+K(d7,d4)+K(d8,d4)+K(d7,d5)+K(d8,d5)+K(d7,d6)+K(d8,d6), K(d7,d7)+K(d7,d8)+K(d8,d7)+K(d8,d8)];
    U(d1,:) = ufixed;
    Fr = -[K(d2,d1); K(d3,d1)+K(d4,d1)+K(d5,d1)+K(d6,d1); K(d7,d1)+K(d8,d1)]*U(d1,:)-[K(d2,d4);K(d3,d4)+K(d4,d4)+K(d5,d4)+K(d6,d4);K(d7,d4)+K(d8,d4)]*w1-[K(d2,d5);K(d3,d5)+K(d4,d5)+K(d5,d5)+K(d6,d5);K(d7,d5)+K(d8,d5)]*w2-[K(d2,d6);K(d3,d6)+K(d4,d6)+K(d5,d6)+K(d6,d6);K(d7,d6)+K(d8,d6)]*w3-[K(d2,d8);K(d3,d8)+K(d4,d8)+K(d5,d8)+K(d6,d8);K(d7,d8)+K(d8,d8)]*w4;
    
    U([d2,d3,d7],:) = Kr\Fr;
    
    U(d4,:) = U(d3,:)+w1;
    U(d5,:) = U(d3,:)+w2;
    U(d6,:) = U(d3,:)+w3;
    U(d8,:) = U(d7,:)+w4;
    
    % PLOT DEFORMATIONS
    
%     for i=1:6
%         figure()
%         plot_def(U(:,i),nelx,nely,nelz);
%     end 

    
       
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
 %figure
%display_3D(xPhys);

% [PSNR,MSE,MAXERR,L2RAT] = measerr(X,XAPP)