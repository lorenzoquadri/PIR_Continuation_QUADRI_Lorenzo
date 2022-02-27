% OPTIMISED CODE
% [tens,obj,micro]=unitCell8tz_3D_optimised(9,9,9,0.5,3,1.5,1,0,0,0,0.2,200,0,0,0,0,0,1,0.5);
%% PERIODIC MATERIAL MICROSTRUCTURE DESIGN
function [tens,obj,micro]=unitCell8tz_3D_optimised(nelx,nely,nelz,density,penal,rmin,ft,ftBC,eta,beta,move,maxit, angle1,angle2,angle3,cubicity21,cubicity31,initDes,transmiLim)
tic
%density : 0 for void, 1 for full material
%angle : 0 for 0 rad, 1 for pi/4 rads
%cubicity : 0 for only one privileged direction, 1 for cubic material
%initDes : initial design : 1=all 0s; 2= all volfrac; 3= all 1; ...
%transmiLim : threshold to be considered as a transmission point
volfrac=density;
cubicity21=sqrt(cubicity21);
cubicity31=sqrt(cubicity31);
angle1=angle1*pi/4;
angle2=angle2*pi/4;
angle3=angle3*pi/4;

addpath('stenglibmaster\Fast')


%% MATERIAL PROPERTIES AND CONTINUATION PARAMETERS
E0 = 1;                                                                    % Young modulus of solid
Emin = 1e-9;                                                               % Young modulus of "void"
nu = 0.3;                                                                  % Poisson ratio
% %
% penalCnt = { 1, 1, 25, 0.25 };                                             % continuation scheme on penal
% betaCnt  = { 1, 1, 25,    2 };                                             % continuation scheme on beta
if ftBC == 'N', bcF = 'symmetric'; else, bcF = 0; end                      % filter BC selector
% %
%% PREPARE FINITE ELEMENT ANALYSIS
nEl = nelx * nely * nelz;                                                  % number of elements          #3D#
nodenrs = int32( reshape( 1 : ( 1 + nelx ) * ( 1 + nely ) * ( 1 + nelz ), ...
    1 + nely, 1 + nelx, 1 + nelz ) );                                      % nodes numbering             #3D#
cVec = reshape( 3 * nodenrs( 1 : nely, 1 : nelx, 1 : nelz ) + 1, nEl, 1 ); %                             #3D#
cMat = cVec+int32( [0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1 3*(nelx+1)*(nely+1)+[0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1]]);             % connectivity matrix         #3D#
nDof = ( 1 + nely ) * ( 1 + nelz ) * ( 1 + nelx ) * 3;                     % total number of DOFs        #3D#
[ sI, sII ] = deal( [ ] );
for j = 1 : 24
    sI = cat( 2, sI, j : 24 );
    sII = cat( 2, sII, repmat( j, 1, 24 - j + 1 ) );
end
[ iK , jK ] = deal( cMat( :,  sI )', cMat( :, sII )' );
Iar = sort( [ iK( : ), jK( : ) ], 2, 'descend' ); clear iK jK              % reduced assembly indexing
Ke = keGen_optimised(nu);                                                            % elemental stiffness matrix
Ke0( tril( ones( 24 ) ) == 1 ) = Ke';
Ke0 = reshape( Ke0, 24, 24 );
Ke0 = Ke0 + Ke0' - diag( diag( Ke0 ) );                                    % recover full matrix

% %
% %% DEFINE IMPLICIT FUNCTIONS
% prj = @(v,eta,beta) (tanh(beta*eta)+tanh(beta*(v(:)-eta)))./...
%     (tanh(beta*eta)+tanh(beta*(1-eta)));                                   % projection
% deta = @(v,eta,beta) - beta * csch( beta ) .* sech( beta * ( v( : ) - eta ) ).^2 .* ...
%     sinh( v( : ) * beta ) .* sinh( ( 1 - v( : ) ) * beta );                % projection eta-derivative
% dprj = @(v,eta,beta) beta*(1-tanh(beta*(v-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));% proj. x-derivative
% cnt = @(v,vCnt,l) v+(l>=vCnt{1}).*(v<vCnt{2}).*(mod(l,vCnt{3})==0).*vCnt{4};
% %
%% PREPARE FILTER

% [dy,dx,dz]=meshgrid(-ceil(rmin)+1:ceil(rmin)-1,...
%     -ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1 );
% h = max( 0, rmin - sqrt( dx.^2 + dy.^2 + dz.^2 ) );                        % conv. kernel #3D#
% iH = ones(nelx*nely*nelz*(2*(ceil(rmin)-1)+1)^2,1);
% jH = ones(size(iH));
% sH = zeros(size(iH));
% % H = sparse(iH,jH,sH);
% Hs = imfilter( ones( nely, nelx, nelz ), h, bcF );                         % matrix of weights (filter)  #3D#
% dHs = Hs;

% %
% %% ALLOCATE AND INITIALIZE OTHER PARAMETERS
% [ dsK, dV ] = deal( zeros( nEl, 1 ) );                                  % initialize vectors
% dV( :, 1 ) = 1/nEl/volfrac;                                              % derivative of volume
% [xOld, change, loop, U ] = deal(1, 1, 0, zeros(3*nEl,6) );       % old x, x change, it. counter, U
% %
% 
iH = ones(nelx*nely*nelz*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
dHs = Hs;
%% PERIODIC BOUNDARY CONDITIONS
e0 = eye(6);
ufixed = zeros(24,6);
U = zeros(nDof,6);
alldofs = int32(1:nDof);
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
x=initDesMore8tz_3D(nelx,nely,nelz,volfrac,initDes,1);
xPhys = x;
change = 1;
loop = 0;
inLoop=1;
%%START ITERATION
while (change > 1e-6 && loop < maxit && inLoop==1) || inLoop==2
    xPhys=(xPhys+flip(fliplr(flipud(xPhys)),3))/2; %make design central symmetric
    loop = loop+1;
    
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
  
    w1 = [repmat(ufixed(10:12,:),(nfarthestUsedy0x0z0-1)+(nelx+1-nfarthestUsedy0xMz0),1);repmat(ufixed(4:6,:),(nfarthestUsedx0y0z0-1)+(nely+1-nfarthestUsedx0yMz0),1);  repmat(ufixed(4:6,:),(nfarthestUsedz0x0y0-1)+(nelx+1-nfarthestUsedzMx0y0),1)];
    w2 = [repmat(ufixed(22:24,:),(nfarthestUsedy0x0z0-1)+(nelx+1-nfarthestUsedy0xMz0),1);repmat(ufixed(16:18,:),(nfarthestUsedx0y0z0-1)+(nely+1-nfarthestUsedx0yMz0),1);  repmat(ufixed(7:9,:),(nfarthestUsedz0x0y0-1)+(nelx+1-nfarthestUsedzMx0y0),1)];
    w3 = [repmat(ufixed(13:15,:),(nfarthestUsedy0x0z0-1)+(nelx+1-nfarthestUsedy0xMz0),1);repmat(ufixed(13:15,:),(nfarthestUsedx0y0z0-1)+(nely+1-nfarthestUsedx0yMz0),1);  repmat(ufixed(10:12,:),(nfarthestUsedz0x0y0-1)+(nelx+1-nfarthestUsedzMx0y0),1)];
    w4 = [repmat(ufixed(4:6,:),(nfarthestUsedx0y0z0-1+nely+1-nfarthestUsedx0yMz0)*(nfarthestUsedz0x0y0-1+nelz+1-nfarthestUsedzMx0y0),1); repmat(ufixed(10:12,:),(nfarthestUsedy0x0z0-1+nelx+1-nfarthestUsedy0xMz0)*(nfarthestUsedz0x0y0-1+nely+1-nfarthestUsedzMx0y0),1); repmat(ufixed(13:15,:),(nfarthestUsedx0y0z0-1+nely+1-nfarthestUsedx0yMz0)*(nfarthestUsedy0x0z0-1+nelx+1-nfarthestUsedy0xMz0),1)];
   
    %% FE-ANALYSIS
    sK = ( Emin + xPhys(:).^penal * ( E0 - Emin ) );
    sK = reshape( Ke( : ) * sK', length( Ke ) * nEl, 1 );                   % transform sK into matrix
    K = fsparse( Iar( :, 1 ), Iar( :, 2 ), sK, [ nDof, nDof ] );  
    K = K + K' - diag( diag( K ) ); K = (K+K')/2;                                   % recover full matrix
    
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
        if loop==1
            for i=1:6
                figure()
                plot_def(U(:,i),nelx,nely,nelz);
            end
        end
    
     
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    for i = 1:6
        for j = 1:6
            U1 = U(:,i); U2 = U(:,j);
            qe{i,j} = reshape(sum((U1(cMat)*Ke0).*U2(cMat),2),nely,nelx,nelz)/(nelx*nely*nelz);
            Q(i,j) = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*qe{i,j})));
            dQ{i,j} = penal*(E0-Emin)*xPhys.^(penal-1).*qe{i,j};
        end
    end
    Q2=rotateTensorMatrix_3D(Q,angle1,angle2,angle3);
    c = -(1-cubicity21*(0.5-1/6*cubicity31)-cubicity31*(0.5-1/6*cubicity21))*Q2(1,1)-cubicity21*(0.5-1/6*cubicity31)*Q2(2,2)-cubicity31*(0.5-1/6*cubicity21)*Q2(3,3);
    dQ2=rotateTensorCells_3D(dQ,angle1,angle2,angle3);
    dc = -(1-cubicity21*(0.5-1/6*cubicity31)-cubicity31*(0.5-1/6*cubicity21))*dQ2{1,1}-cubicity21*(0.5-1/6*cubicity31)*dQ2{2,2}-cubicity31*(0.5-1/6*cubicity21)*dQ2{3,3};
    dv = ones(nely,nelx,nelz);
    %% FILTERING/MODIFICATION OF SENSITIVITIES
%     %
%     xPhys = imfilter( reshape( x, nely, nelz, nelx ), h, bcF ) ./ Hs;       % filtered field              #3D#
%                                              
%     if ft > 1                              % compute optimal eta* with Newton
%         f = ( mean( prj( xPhys(:), eta, beta ) ) - volfrac )  * (ft == 3);      % function (volume)
%         while abs( f ) > 1e-6           % Newton process for finding opt. eta
%             eta = eta - f / mean( deta( xPhys(:), eta, beta ) );
%             f = mean( prj( xPhys(:), eta, beta ) ) - volfrac;
%         end
%         dHs = Hs ./ reshape( dprj( xPhys(:), eta, beta ), nely, nelz, nelx );   % sensitivity modification    #3D#
%         xPhys(:) = prj( xPhys(:), eta, beta );                                     % projected (physical) field
%     end
%     change = norm( xPhys(:) - xOld(:) ) ./ nEl;
%     xOld = xPhys; 
%     %
    if ft==1 %sensitivity
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
        %dc = imfilter( reshape( dc, nely, nelz, nelx ) ./ dHs, h, bcF );         % filter objective sens.      #3D#
    elseif ft == 2 % density
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
%     xT = x( : );
%     [ xU, xL ] = deal( xT + move, xT - move );                               % current upper and lower bound
%     ocP = xT .* sqrt( - dc( : ) ./ dV0( : ) );                           % constant part in resizing rule
%     l = [ 0, mean( ocP ) / volfrac ];                                        % initial estimate for LM
%     while ( l( 2 ) - l( 1 ) ) / ( l( 2 ) + l( 1 ) ) > 1e-4                   % OC resizing rule
%         lmid = 0.5 * ( l( 1 ) + l( 2 ) );
%         x( : ) = max( max( min( min( ocP / lmid, xU ), 1 ), xL ), 0 );
%         if mean( x ) > volfrac, l( 1 ) = lmid; else, l( 2 ) = lmid; end
%     end
%     [penal,beta] = deal(cnt(penal,penalCnt,loop), cnt(beta,betaCnt,loop));   % apply conitnuation on parameters
%     
%     if inLoop==2
%         inLoop=0;
%         obj=c;
%         tens=Q;
%     elseif (change <= 1e-6 || loop >=200) && inLoop==1
%         inLoop=2;
%         micro=xPhys;
%     end

    l1 = 0; l2 = 1e9; move = 0.2;
    while(l2-l1 > 1e-9)
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(max(0,-dc./dv/lmid))))));
        %xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        if ~isreal(xnew)
            debug=1
        end
        if ft ==1
            xPhys = xnew;
        elseif ft == 2
            xPhys(:) = (H*xnew(:))./Hs;
        end
        if mean(xPhys(:)) > volfrac, l1 = lmid; else l2 = lmid; end
    end
    if ~isreal(xPhys)
        debug=1;
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    
    if inLoop==2
        inLoop=0;
        obj=c;
        tens=Q;
    elseif (change <= 0.01 || loop >=200) && inLoop==1
        inLoop=2;
        micro=xPhys;
    end
    
    
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xPhys(:)),change);
    clf;
    %display_3D(xPhys);
    toc
end

