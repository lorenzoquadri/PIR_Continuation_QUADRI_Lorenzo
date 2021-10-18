% [cO,xdensO,xcos1O,xsin1O,xcos2O,xsin2O,xcos3O,xsin3O,xcub21O,xcub31O]=topMulti_3D(10,5,1,0.5,"volfrac","Canti")
function [cO,xdensO,xcos1O,xsin1O,xcos2O,xsin2O,xcos3O,xsin3O,xcub21O,xcub31O]=topMulti_3D(nelx,nely,nelz,volfrac,initialDesign,problem)
% USER-DEFINED MODEL PARAMETERS
%nelx : number of cells in horizontal direction
%nely : number of cells in vertical direction
%nelz: number of cells in depth direction
%volfrac : global volume fraction
rmin = 1.5; %filtering radius
fsum=1.0; %force value
xMin = 0; %minimum cell density
xMax = 1; %maximum cell density


global B database sig;
sig=0.15; %gaussion kernel radius
B = func_B();
load('database4-4-4-4-4-4.mat'); % cell elastic tensor database
database=dbMat;



% USER-DEFINED LOOP PARAMETERS
maxloopaftermin=5; % Maximum number of iterations without a new global minimum
maxloop=100; % Maximum number of iterations
tolx = 0.001; % Terminarion criterion

switch problem
    %     case 'MBB'
    %     % USER-DEFINED LOAD DOFs
    %     [il1,jl1,kl1] = meshgrid(0,nely,0:nelz);                  % Coordinates
    %     loadnid = kl1*(nelx+1)*(nely+1)+il1*(nely+1)+(nely+1-jl1); % Node IDs
    %     loaddof = 2*loadnid(:) ; % DOFs
    %     % USER-DEFINED SUPPORT FIXED DOFs
    %     [if1,jf1,kf1] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
    %     fixednid_1 = kf1*(nelx+1)*(nely+1)+if1*(nely+1)+(nely+1-jf1); % Node IDs
    %     [if2,jf2,kf2] = meshgrid(nelx,0,0:nelz);                  % Coordinates
    %     fixednid_2 = kf2*(nelx+1)*(nely+1)+if2*(nely+1)+(nely+1-jf2); % Node IDs
    %     fixeddof = [3*fixednid_1(:)-2; 3*fixednid_1(:); 3*fixednid_2(:)-1; 3*fixednid_2(:)]; % DOFs
    %     % USER-DEFINED ACTIVE ELEMENTS
    %     activeelts=ones(nelx*nely,1); %%%
    case 'Canti'
        % USER-DEFINED LOAD DOFs
        [il1,jl1,kl1] = meshgrid(nelx,0,0:nelz);                  % Coordinates
        loadnid = kl1*(nelx+1)*(nely+1)+il1*(nely+1)+(nely+1-jl1); % Node IDs
        loaddof = 2*loadnid(:) ; % DOFs
        % USER-DEFINED SUPPORT FIXED DOFs
        [if1,jf1,kf1] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
        fixednid_1 = kf1*(nelx+1)*(nely+1)+if1*(nely+1)+(nely+1-jf1); % Node IDs
        fixeddof = [3*fixednid_1(:)-2; 3*fixednid_1(:)-1; 3*fixednid_1(:)]; % DOFs
        % USER-DEFINED ACTIVE ELEMENTS
        activeelts=ones(nelx*nely*nelz,1); %%%
        %     case 'Lshape'
        %     % USER-DEFINED LOAD DOFs
        %     [il1,jl1,kl1] = meshgrid(nelx,0,0:nelz);                  % Coordinates
        %     loadnid = kl1*(nelx+1)*(nely+1)+il1*(nely+1)+(nely+1-jl1); % Node IDs
        %     loaddof = 2*loadnid(:) ; % DOFs
        %     % USER-DEFINED SUPPORT FIXED DOFs
        %     [if1,jf1,kf1] = meshgrid(0:fix(nelx/2),nely,0:nelz);                  % Coordinates
        %     fixednid_1 = kf1*(nelx+1)*(nely+1)+if1*(nely+1)+(nely+1-jf1); % Node IDs
        %     fixeddof = [3*fixednid_1(:)-2; 3*fixednid_1(:)-1; 3*fixednid_1(:)]; % DOFs
        %     % USER-DEFINED ACTIVE ELEMENTS
        %     emptyelts=(nelx/2)*(nely)+1:(nelx)*(nely);
        %     emptyelts=reshape(emptyelts, nely,nelx/2);
        %     emptyelts=emptyelts(1:nely/2,:);
        %     emptyelts=emptyelts(:);
        %     activeelts=ones(nelx*nely,1);
        %     activeelts(emptyelts)=0;
end

% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-fsum,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
volfrac=volfrac*mean(activeelts);

nodenrs = reshape(1:(1+nelx)*(1+nely)*(1+nelz),1+nely,1+nelx,1+nelz);
edofVec = reshape(3*nodenrs(1:end-1,1:end-1,1:end-1)+1, nelx*nely*nelz,1);
edofMat = repmat(edofVec,1,24)+repmat([0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1 3*(nelx+1)*(nely+1)+[0 1 2 3*nely+[3 4 5 0 1 2] -3 -2 -1]],nelx*nely*nelz,1);
iK = reshape(kron(edofMat,ones(24,1))',576*nelx*nely*nelz,1);
jK = reshape(kron(edofMat,ones(1,24))',576*nelx*nely*nelz,1);

%% PREPARE FILTER
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


%% INITIALIZE ITERATION
if initialDesign=="top88"
    switch problem
        %         case 'MBB'
        %         xdens = top88DesignMBB(nelx,nely,volfrac,2,1.2,2); xcos = repmat(1, [nely, nelx]); xsin = repmat(1, [nely, nelx]); xcub = repmat(0.5, [nely, nelx]);
        case 'Canti'
            xdens = top88DesignCanti(nelx,nely,nelz,volfrac,2,1.2,2); xcos1 = ones(nely,nelx,nelz); xsin1 = ones(nely,nelx,nelz); xcos2 = ones(nely,nelx,nelz); xsin2 = ones(nely,nelx,nelz); xcos3 = ones(nely,nelx,nelz); xsin3 = ones(nely,nelx,nelz); xcub21 = 0.5*ones(nely,nelx,nelz); xcub31 = 0.5*ones(nely,nelx,nelz);
            %         case 'Lshape'
            %         xdens = top88DesignL(nelx,nely,volfrac/mean(activeelts),2,1.2,2); xcos = repmat(1, [nely, nelx]); xsin = repmat(1, [nely, nelx]); xcub = repmat(0.5, [nely, nelx]);
    end
elseif initialDesign=="volfrac"
    xdens = repmat(volfrac, [nely, nelx, nelz]); xcos1 = ones(nely,nelx,nelz); xsin1 = ones(nely,nelx,nelz); xcos2 = ones(nely,nelx,nelz); xsin2 = ones(nely,nelx,nelz); xcos3 = ones(nely,nelx,nelz); xsin3 = ones(nely,nelx,nelz); xcub21 = 0.5*ones(nely,nelx,nelz); xcub31 = 0.5*ones(nely,nelx,nelz);
end
xdensPhys = xdens; xcos1Phys = xcos1; xsin1Phys = xsin1; xcos2Phys = xcos2; xsin2Phys = xsin2; xcos3Phys = xcos3; xsin3Phys = xsin3; xcub21Phys = xcub21; xcub31Phys = xcub31;
loop = 0; change = 1; loopaftermin=0;

% INITIALIZE MMA OPTIMIZER
m = 1; n = 9*nele;
xmin = [xMin*ones(nele,1); zeros(nele,1); zeros(nele,1); zeros(nele,1); zeros(nele,1); zeros(nele,1); zeros(nele,1); zeros(nele,1); zeros(nele,1)]; % Column vector with the lower bounds for the macro-variables.
xmax = [xMax*ones(nele,1); ones(nele,1); ones(nele,1); ones(nele,1); ones(nele,1); ones(nele,1); ones(nele,1); ones(nele,1); ones(nele,1)]; % Column vector with the upper bounds for the macro-variables.
xval = [xdensPhys(:); xcos1Phys(:); xsin1Phys(:); xcos2Phys(:); xsin2Phys(:); xcos3Phys(:); xsin3Phys(:); xcub21Phys(:); xcub31Phys(:)]; % macro-variables
xold1 = xval(:); xold2 = xold1(:);
low = ones(n,1); upp = ones(n,1);
a0 = 1; a_mma = zeros(m,1); c_mma = 5000*ones(m,1); d_mma = zeros(m,1);

%INITIALIZE GLOBAL OPTIMUM
xdensO=zeros(nely,nelx,nelz);
xcos1O=zeros(nely,nelx,nelz);
xsin1O=zeros(nely,nelx,nelz);
xcos2O=zeros(nely,nelx,nelz);
xsin2O=zeros(nely,nelx,nelz);
xcos3O=zeros(nely,nelx,nelz);
xsin3O=zeros(nely,nelx,nelz);
xcub21O=zeros(nely,nelx,nelz);
xcub31O=zeros(nely,nelx,nelz);
cO=inf;
ceO=zeros(nely,nelx,nelz);

% START ITERATION
while change > tolx && loop < maxloop && loopaftermin < maxloopaftermin
    loop = loop+1;
    loopaftermin = loopaftermin+1;
    % FE-ANALYSIS AND SENSITIVITY ANALYSIS
    [K_cell, K_dxdens_cell, K_dxcos1_cell, K_dxsin1_cell, K_dxcos2_cell, K_dxsin2_cell,K_dxcos3_cell, K_dxsin3_cell, K_dxcub21_cell, K_dxcub31_cell] = arrayfun(@KE_matrix, xdensPhys(:)', xcos1Phys(:)', xsin1Phys(:)', xcos2Phys(:)', xsin2Phys(:)', xcos3Phys(:)', xsin3Phys(:)', xcub21Phys(:)',xcub31Phys(:)', 'un', 0);
    
    KALL = reshape(1/8*cell2mat(K_cell), [24*24, nele]);
    sK = reshape(KALL, 24*24*nele, 1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    if max(max(abs(U)))>2000
        U=2000*U/max(max(abs(U))); % rescale U if it is too big for MMA to handle properly
    end
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    F_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_cell,'un', 0)'); %1/8
    ce = reshape(sum(F_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    c = sum(sum(sum(ce)));
    %SAVE GLOBAL OPTIMUM
    if mean(xdensPhys(:)) <= volfrac && c < cO
        xdensO=xdens; % add other terms
        xcos1O=xcos1;
        xsin1O=xsin1;
        xcos2O=xcos2;
        xsin2O=xsin2;
        xcos3O=xcos3;
        xsin3O=xsin3;
        xcub21O=xcub21;
        xcub31O=xcub31;
        cO=c;
        ceO=ce;
        loopaftermin = 0;
    end
    
    F_dxdens_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxdens_cell,'un', 0)');
    ce_dxdens = reshape(sum(F_dxdens_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xdens = -ce_dxdens;
    
    F_dxcos1_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxcos1_cell,'un', 0)');
    ce_dxcos1 = reshape(sum(F_dxcos1_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xcos1 = -ce_dxcos1;
    
    F_dxsin1_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxsin1_cell,'un', 0)');
    ce_dxsin1 = reshape(sum(F_dxsin1_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xsin1 = -ce_dxsin1;
    
    F_dxcos2_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxcos2_cell,'un', 0)');
    ce_dxcos2 = reshape(sum(F_dxcos2_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xcos2 = -ce_dxcos2;
    
    F_dxsin2_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxsin2_cell,'un', 0)');
    ce_dxsin2 = reshape(sum(F_dxsin2_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xsin2 = -ce_dxsin2;
    
    F_dxcos3_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxcos3_cell,'un', 0)');
    ce_dxcos3 = reshape(sum(F_dxcos3_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xcos3 = -ce_dxcos3;
    
    F_dxsin3_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxsin3_cell,'un', 0)');
    ce_dxsin3 = reshape(sum(F_dxsin3_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xsin3 = -ce_dxsin3;
    
    F_dxcub21_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxcub21_cell,'un', 0)');
    ce_dxcub21 = reshape(sum(F_dxcub21_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xcub21 = -ce_dxcub21;
    
    F_dxcub31_nodes = cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dxcub31_cell,'un', 0)');
    ce_dxcub31 = reshape(sum(F_dxcub31_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
    dc_xcub31 = -ce_dxcub31;
    
    dv_x = ones(nely,nelx,nelz);
    
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    dc_xdens(:) = H*(xdens(:).*dc_xdens(:))./Hs./max(1e-4,xdens(:));
    
    % MMA OPTIMIZATION METHOD
    %
    f0val = c; df0dx = [dc_xdens(:).*activeelts; dc_xcos1(:).*activeelts; dc_xsin1(:).*activeelts; dc_xcos2(:).*activeelts; dc_xsin2(:).*activeelts; dc_xcos3(:).*activeelts; dc_xsin3(:).*activeelts; dc_xcub21(:).*activeelts; dc_xcub31(:).*activeelts];
    fval = sum(xdensPhys(:))/(volfrac*nele) - 1;
    dfdx = [(dv_x(:).*activeelts)'/(volfrac*nele), zeros(1,nele), zeros(1,nele), zeros(1,nele), zeros(1,nele), zeros(1,nele), zeros(1,nele), zeros(1,nele), zeros(1,nele)];
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a_mma,c_mma,d_mma);
    xold2 = xold1; xold1 = xval; change = max(abs(xmma-xval)); xval = xmma;
    
    xdensnew = reshape(xval(1:nele), nely, nelx, nelz);
    xcos1new = reshape(xval(nele+1:2*nele), nely, nelx, nelz);
    xsin1new = reshape(xval(2*nele+1:3*nele), nely, nelx, nelz);
    xcos2new = reshape(xval(3*nele+1:4*nele), nely, nelx, nelz);
    xsin2new = reshape(xval(4*nele+1:5*nele), nely, nelx, nelz);
    xcos3new = reshape(xval(5*nele+1:6*nele), nely, nelx, nelz);
    xsin3new = reshape(xval(6*nele+1:7*nele), nely, nelx, nelz);
    xcub21new = reshape(xval(7*nele+1:8*nele), nely, nelx, nelz);
    xcub31new = reshape(xval(8*nele+1:9*nele), nely, nelx, nelz);
    
    % FILTERING AND MODIFICATION OF VARIABLES
    xdensnew(:) = xdensnew(:).*activeelts;
    xcos1new(:) = (H*xcos1new(:))./Hs; xcos1new(xcos1new > 1.0) = 1.0;
    xsin1new(:) = (H*xsin1new(:))./Hs; xsin1new(xsin1new > 1.0) = 1.0;
    xcos2new(:) = (H*xcos2new(:))./Hs; xcos2new(xcos2new > 1.0) = 1.0;
    xsin2new(:) = (H*xsin2new(:))./Hs; xsin2new(xsin2new > 1.0) = 1.0;
    xcos3new(:) = (H*xcos3new(:))./Hs; xcos3new(xcos3new > 1.0) = 1.0;
    xsin3new(:) = (H*xsin3new(:))./Hs; xsin3new(xsin3new > 1.0) = 1.0;
    xcub21new(:) = H*(xcub21new(:)./Hs); xcub21new(xcub21new > 1.0) = 1.0;
    xcub31new(:) = H*(xcub31new(:)./Hs); xcub31new(xcub31new > 1.0) = 1.0;
    xdens = xdensnew; xcos1 = xcos1new; xsin1 = xsin1new; xcos2 = xcos2new; xsin2 = xsin2new; xcos3 = xcos3new; xsin3 = xsin3new; xcub21 = xcub21new; xcub31 = xcub31new;
    
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f \n',loop,c,mean(xdensPhys(:)),change);
    xdensPhys = xdens; xcos1Phys = xcos1; xsin1Phys = xsin1; xcos2Phys = xcos2; xsin2Phys = xsin2; xcos3Phys = xcos3; xsin3Phys = xsin3; xcub21Phys = xcub21; xcub31Phys = xcub31;
    
    figure(1)
    title('Density of cells')
    display_3D(xdensPhys);
    
    figure(2)
    title('cos1 of cells')
    display_3D(xcos1Phys.*reshape(activeelts,nely,nelx,nelz));
    figure(3)
    title('sin1 of cells')
    display_3D(xsin1Phys.*reshape(activeelts,nely,nelx,nelz));
    
    figure(4)
    title('cos2 of cells')
    display_3D(xcos2Phys.*reshape(activeelts,nely,nelx,nelz));
    figure(5)
    title('sin2 of cells')
    display_3D(xsin2Phys.*reshape(activeelts,nely,nelx,nelz));
    
    figure(6)
    title('cos3 of cells')
    display_3D((1-xcos3Phys).*reshape(activeelts,nely,nelx,nelz));
    figure(7)
    title('sin3 of cells')
    display_3D(xsin3Phys.*reshape(activeelts,nely,nelx,nelz));
    
    figure(8)
    title('cub21 of cells')
    display_3D((1-xcub21Phys).*reshape(activeelts,nely,nelx,nelz));
    figure(9)
    title('cub31 of cells')
    display_3D(xcub31Phys.*reshape(activeelts,nely,nelx,nelz));end
end

%% GEOMETRIC MATRIX B USED IN STIFFNESS MATRIX CALCULATION
function B = func_B()
syms s t n;
x1 = 0; y1 = 0; z1 = 0; x2 = 1; y2 = 0; z2 = 0;
x3 = 1; y3 = 1; z3 = 0; x4 = 0; y4 = 1; z4 = 0;
x5 = 0; y5 = 0; z5 = 1; x6 = 1; y6 = 0; z6 = 1;
x7 = 1; y7 = 1; z7 = 1; x8 = 0; y8 = 1; z8 = 1;
N1 = (1+s)*(1-t)*(1-n)/8; N2 = (1+s)*(1+t)*(1-n)/8;
N3 = (1-s)*(1+t)*(1-n)/8; N4 = (1-s)*(1-t)*(1-n)/8;
N5 = (1+s)*(1-t)*(1+n)/8; N6 = (1+s)*(1+t)*(1+n)/8;
N7 = (1-s)*(1+t)*(1+n)/8; N8 = (1-s)*(1-t)*(1+n)/8;
x = N1*x1 + N2*x2 + N3*x3 + N4*x4 + N5*x5 + N6*x6 + N7*x7 + N8*x8;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4 + N5*y5 + N6*y6 + N7*y7 + N8*y8;
z = N1*z1 + N2*z2 + N3*z3 + N4*z4 + N5*z5 + N6*z6 + N7*z7 + N8*z8;
J = [diff(x,s), diff(y,s), diff(z,s);
    diff(x,t), diff(y,t), diff(z,t);
    diff(x,n), diff(y,n), diff(z,n)];
Jdet = det(J);
a = diff(y,t)*diff(z,n)-diff(z,t)*diff(y,n);
b = diff(y,s)*diff(z,n)-diff(z,s)*diff(y,n);
c = diff(y,s)*diff(z,t)-diff(z,s)*diff(y,t);
d = diff(x,t)*diff(z,n)-diff(z,t)*diff(x,n);
e = diff(x,s)*diff(z,n)-diff(z,s)*diff(x,n);
f = diff(x,s)*diff(z,t)-diff(z,s)*diff(x,t);
g = diff(x,t)*diff(y,n)-diff(y,t)*diff(x,n);
h = diff(x,s)*diff(y,n)-diff(y,s)*diff(x,n);
l = diff(x,s)*diff(y,t)-diff(y,s)*diff(x,t);
Ns = [N1, N2, N3, N4, N5, N6, N7, N8];
Bs = sym(zeros(6,3,8));
for i = 1:8
    Bs(:,:,i) = [ a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i),n), 0, 0; %ho messo ultima n
        0, -d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i),n), 0; %%ho messo ultima n
        0, 0, g*diff(Ns(i),t)-h*diff(Ns(i), t)+l*diff(Ns(i),n);
        -d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i),n), a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i),n), 0;
        0, g*diff(Ns(i),s)-h*diff(Ns(i),t)+l*diff(Ns(i),n), -d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i),n);
        g*diff(Ns(i),s)-h*diff(Ns(i),t)+l*diff(Ns(i),n), 0, a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i),n)]/Jdet;
end
B = [Bs(:,:,1),Bs(:,:,2),Bs(:,:,3),Bs(:,:,4),Bs(:,:,5),Bs(:,:,6),Bs(:,:,7),Bs(:,:,8)];
B = matlabFunction(B);
end
% CALCULATION OF STIFFNESS MATRIX AND ITS FOUR PARTIAL DERIVATIVES
% The stiffness matrix and derivatives depends on the cell microstructure, defined by the
% macro-variables
function [KE, KE_dxdens, KE_dxcos1, KE_dxsin1, KE_dxcos2, KE_dxsin2,KE_dxcos3, KE_dxsin3, KE_dxcub21, KE_dxcub31] = KE_matrix(xdens, xcos1, xsin1, xcos2, xsin2, xcos3, xsin3, xcub21, xcub31)
global B database sig;
%% get point's tensor and derivatives
% density
xdensd=xdens+0.01; %xdens+delta used to get the partial derivative approximation
negdifdens=1;
if xdensd>1 %if right partial derivative isn't accessible, get left partial derivative
    xdensd=xdensd-0.02;
    negdifdens=-1;
end

% cubicity
xcub21d=xcub21+0.01; %xcub+delta used to get the partial derivatives
negdifcub21=1;
if xcub21d>1 %if right partial derivative isn't accessible, get left partial derivative
    xcub21d=xcub21d-0.02;
    negdifcub21=-1;
end
xcub31d=xcub31+0.01; %xcub+delta used to get the partial derivatives
negdifcub31=1;
if xcub31d>1 %if right partial derivative isn't accessible, get left partial derivative
    xcub31d=xcub31d-0.02;
    negdifcub31=-1;
end

%derive xo1r from xcos and xsin.
cosalpha1=2*xcos1-1; sinalpha1=2*xsin1-1;
xo1r=atan(sinalpha1/cosalpha1)/pi; %Here, xor is in [-1,1], representing an orientation angle in [-pi,pi]
%put the orientation angle in [0,pi]
if xo1r<0
    xo1r=xo1r+1;
end
xo1rd=xo1r+0.01; %xcor+delta used to get the partial derivatives
negdifor1=1;
if xo1rd>1 %if right partial derivative isn't accessible, get left partial derivative
    xo1rd=xo1rd-0.02;
    negdifor1=-1;
end

%derive xo2r from xcos and xsin.
cosalpha2=2*xcos2-1; sinalpha2=2*xsin2-1;
xo2r=atan(sinalpha2/cosalpha2)/pi; %Here, xor is in [-1,1], representing an orientation angle in [-pi,pi]
%put the orientation angle in [0,pi]
if xo2r<0
    xo2r=xo2r+1;
end
xo2rd=xo2r+0.01; %xcor+delta used to get the partial derivatives
negdifor2=1;
if xo2rd>1 %if right partial derivative isn't accessible, get left partial derivative
    xo2rd=xo2rd-0.02;
    negdifor2=-1;
end

%derive xo3r from xcos and xsin.
cosalpha3=2*xcos3-1; sinalpha3=2*xsin3-1;
xo3r=atan(sinalpha3/cosalpha3)/pi; %Here, xor is in [-1,1], representing an orientation angle in [-pi,pi]
%put the orientation angle in [0,pi]
if xo3r<0
    xo3r=xo3r+1;
end
xo3rd=xo3r+0.01; %xcor+delta used to get the partial derivatives
negdifor3=1;
if xo3rd>1 %if right partial derivative isn't accessible, get left partial derivative
    xo3rd=xo3rd-0.02;
    negdifor3=-1;
end

%macrovariables for xi, xi', xi'' and xi'''
xdensv=[xdens, xdensd, xdens, xdens, xdens, xdens, xdens];
xo1rv=[xo1r, xo1r, xo1rd, xo1r, xo1r, xo1r, xo1r];
xo2rv=[xo2r, xo2r, xo2r, xo2rd, xo2r, xo2r, xo2r];
xo3rv=[xo3r, xo3r, xo3r, xo3r, xo3rd, xo3r, xo3r];
xcub21v=[xcub21, xcub21, xcub21, xcub21, xcub21, xcub21d, xcub21];
xcub31v=[xcub31, xcub31, xcub31, xcub31, xcub31, xcub31, xcub31d];

d11=[];
d12=[];
d13=[];
d14=[];
d15=[];
d16=[];
d22=[];
d23=[];
d24=[];
d25=[];
d26=[];
d33=[];
d34=[];
d35=[];
d36=[];
d44=[];
d45=[];
d46=[];
d55=[];
d56=[];
d66=[];

% GET ELASTICITY TENSORS FROM DATABASE METAMODEL FOR Xi, Xi', Xi'' and Xi'''
for i = 1:7
    % find points in matrix within 3 kernel radii
    xdenslim1=max(1,round((xdensv(i)-3*sig)*(size(database,7)-1)+1));
    xdenslim2=min(size(database,7),round((xdensv(i)+3*sig)*(size(database,7)-1)+1));
    xo1rlim1=max(1,round((xo1rv(i)-3*sig)*(size(database,6)-1)+1));
    xo1rlim2=min(size(database,6),round((xo1rv(i)+3*sig)*(size(database,6)-1)+1));
    xo2rlim1=max(1,round((xo2rv(i)-3*sig)*(size(database,5)-1)+1));
    xo2rlim2=min(size(database,5),round((xo2rv(i)+3*sig)*(size(database,5)-1)+1));
    xo3rlim1=max(1,round((xo3rv(i)-3*sig)*(size(database,4)-1)+1));
    xo3rlim2=min(size(database,4),round((xo3rv(i)+3*sig)*(size(database,4)-1)+1));
    xcub21lim1=max(1,round((xcub21v(i)-3*sig)*(size(database,3)-1)+1));
    xcub21lim2=min(size(database,3),round((xcub21v(i)+3*sig)*(size(database,3)-1)+1));
    xcub31lim1=max(1,round((xcub31v(i)-3*sig)*(size(database,2)-1)+1));
    xcub31lim2=min(size(database,2),round((xcub31v(i)+3*sig)*(size(database,2)-1)+1));
    
    Msim=database(:,xcub31lim1:xcub31lim2,xcub21lim1:xcub21lim2,xo3rlim1:xo3rlim2,xo2rlim1:xo2rlim2,xo1rlim1:xo1rlim2,xdenslim1:xdenslim2);
    
    %get distance of each point
    distance=sqrt((Msim(1,:,:,:,:,:,:)-xdensv(i)).^2+(Msim(2,:,:,:,:,:,:)-xo1rv(i)).^2+(Msim(3,:,:,:,:,:,:)-xo2rv(i)).^2+(Msim(4,:,:,:,:,:,:)-xo3rv(i)).^2+(Msim(5,:,:,:,:,:,:)-xcub21v(i)).^2+(Msim(6,:,:,:,:,:,:)-xcub31v(i)).^2);
    %Nadaraya-Watson kernel-weighted average
    gaussFactor = (1/((2*pi*sig^2)^(3/2)))*exp(-((distance).^2)./(2*sig^2));
    gaussFactors = repmat(gaussFactor,21,1,1,1,1,1,1);
    totalGauss=sum(sum(sum(sum(sum(sum(gaussFactor))))));
    pointGauss=sum(sum(sum(sum(sum(sum(Msim(7:end,:,:,:,:,:,:).*gaussFactors,2), 3), 4), 5), 6), 7)/totalGauss; % it's Epred for each i. dij is the vector of the same component o Epred for different i, different variables
    
    
    d11 = [d11,pointGauss(1)];
    d12 = [d12,pointGauss(2)];
    d13 = [d13,pointGauss(3)];
    d14 = [d14,pointGauss(4)];
    d15 = [d15,pointGauss(5)];
    d16 = [d16,pointGauss(6)];
    d22 = [d22,pointGauss(7)];
    d23 = [d23,pointGauss(8)];
    d24 = [d24,pointGauss(9)];
    d25 = [d25,pointGauss(10)];
    d26 = [d26,pointGauss(11)];
    d33 = [d33,pointGauss(12)];
    d34 = [d34,pointGauss(13)];
    d35 = [d35,pointGauss(14)];
    d36 = [d36,pointGauss(15)];
    d44 = [d44,pointGauss(16)];
    d45 = [d45,pointGauss(17)];
    d46 = [d46,pointGauss(18)];
    d55 = [d55,pointGauss(19)];
    d56 = [d56,pointGauss(20)];
    d66 = [d66,pointGauss(21)];
end
%Approximate partial derivatives (21 combinazioni *6)
for i=1:7
    D{i} = [d11(i), d12(i), d13(i), d14(i), d15(i), d16(i),...
        d22(i), d23(i), d24(i), d25(i), d26(i),...
        d33(i), d34(i), d35(i), d36(i),...
        d44(i), d45(i), d46(i),...
        d55(i), d56(i),...
        d66(i)]';
end
D_1=D{1};
negdif=[negdifdens,negdifor1,negdifor2,negdifor3,negdifcub21,negdifcub31];

for i=1:6
    D_dx_array=[];
    D_var=[];
    D_var=D{i+1}; % select with respect to which variable the derivative is calculated
    negdif_var=negdif(i);
    for j=1:21
        D_dx_array(j)=negdif_var*(D_var(j)-D_1(j))/0.01;
    end
    D_dx{i}=D_dx_array;
end

% fill the upper triangular matrix with the derivative elements
loop=1;
D_dxdens_array=D_dx{1};
D_dxor1_array=D_dx{2};
D_dxor2_array=D_dx{3};
D_dxor3_array=D_dx{4};
D_dxcub21_array=D_dx{5};
D_dxcub31_array=D_dx{6};

D_dxdens=array2matrix(D_dxdens_array);
D_dxor1=array2matrix(D_dxor1_array);
D_dxor2=array2matrix(D_dxor2_array);
D_dxor3=array2matrix(D_dxor3_array);
D_dxcub21=array2matrix(D_dxcub21_array);
D_dxcub31=array2matrix(D_dxcub31_array);

%DERIVE STIFFNESS MATRIX FROM ELASTICITY TENSOR
%get gauss points weights and positions
x = [-sqrt(15)/5, 0, sqrt(15)/5];
w = [5/9, 8/9, 5/9];
ke = zeros(24);
ke_dxdens = zeros(24);
ke_dxor1 = zeros(24);
ke_dxor2 = zeros(24);
ke_dxor3 = zeros(24);
ke_dxcub21 = zeros(24);
ke_dxcub31 = zeros(24);

[i, j, k] = meshgrid(1:3, 1:3, 1:3);

for m = 1:27
    ke = ke + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D_dxdens * B(x(i(m)),  x(j(m)), x(k(m)));
    ke_dxor1 = ke_dxor1 + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D_dxor1 * B(x(i(m)), x(j(m)), x(k(m)));
    ke_dxor2 = ke_dxor2 + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D_dxor2 * B(x(i(m)), x(j(m)), x(k(m)));
    ke_dxor3 = ke_dxor3 + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D_dxor3 * B(x(i(m)), x(j(m)), x(k(m)));
    ke_dxcub21 = ke_dxcub21 + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D_dxcub21 * B(x(i(m)), x(j(m)), x(k(m)));
    ke_dxcub31 = ke_dxcub31 + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D_dxcub31 * B(x(i(m)), x(j(m)), x(k(m)));
end

ke_dxor1=ke_dxor1/pi;
ke_dxcosalpha1=ke_dxor1*(-sinalpha1/(cosalpha1^2+sinalpha1^2));
ke_dxcos1=2*ke_dxcosalpha1;
ke_dxsinalpha1=ke_dxor1*(cosalpha1/(cosalpha1^2+sinalpha1^2));
ke_dxsin1=2*ke_dxsinalpha1;

ke_dxor2=ke_dxor2/pi;
ke_dxcosalpha2=ke_dxor2*(-sinalpha2/(cosalpha2^2+sinalpha2^2));
ke_dxcos2=2*ke_dxcosalpha2;
ke_dxsinalpha2=ke_dxor2*(cosalpha2/(cosalpha2^2+sinalpha2^2));
ke_dxsin2=2*ke_dxsinalpha2;

ke_dxor3=ke_dxor3/pi;
ke_dxcosalpha3=ke_dxor3*(-sinalpha3/(cosalpha3^2+sinalpha3^2));
ke_dxcos3=2*ke_dxcosalpha3;
ke_dxsinalpha3=ke_dxor3*(cosalpha3/(cosalpha3^2+sinalpha3^2));
ke_dxsin3=2*ke_dxsinalpha3;

KE = ke; KE_dxdens = ke_dxdens; KE_dxcos1 = ke_dxcos1; KE_dxsin1 = ke_dxsin1; KE_dxcos2 = ke_dxcos2; KE_dxsin2 = ke_dxsin2; KE_dxcos3 = ke_dxcos3; KE_dxsin3 = ke_dxsin3; KE_dxcub21 = ke_dxcub21; KE_dxcub31 = ke_dxcub31;

end