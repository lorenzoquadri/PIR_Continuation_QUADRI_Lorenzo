% A MATLAB CODE FOR CONCURRENT TOPPLOGY OPTIMIZATION BASED ON PILM MODEL BY
% WANG CHUANG AND ZHOU HAN
function PILM()
% USER-DEFINED MODEL PARAMETERS
nelx = 30; nely = 12; nelz = 12;
volfrac = 0.3; rmin = 2.0; fsum = 60;
global B; B = func_B();
% USER-DEFINED LOOP PARAMETERS
maxloop = 200; % Maximum number of iterations
tolx = 0.001; % Terminarion criterion
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(0, nely, 0:nelz);
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1; % DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[if_1,jf_1,kf_1] = meshgrid(0,0:nely,0:nelz); % Coordinates
fixednid_1 = kf_1*(nelx+1)*(nely+1)+if_1*(nely+1)+(nely+1-jf_1); % Node IDs
[if_2,jf_2,kf_2] = meshgrid(nelx*4/5,0,0:nelz); % Coordinates
fixednid_2 = kf_2*(nelx+1)*(nely+1)+if_2*(nely+1)+(nely+1-jf_2);
[if_3,jf_3,kf_3] = meshgrid(nelx*4/5,0,nelz*0.5); % Coordinates
fixednid_3 = kf_3*(nelx+1)*(nely+1)+if_3*(nely+1)+(nely+1-jf_3);
fixeddof = [3*fixednid_1(:)-2; 3*fixednid_2(:)-1; 3*fixednid_3(:)]; % DOFs
% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-fsum/(nelz+1),ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
% Concurrent design of hierarchical structures with three-dimensional parameterized lattice microstructures... 
edofMat = repmat(edofVec,1,24)+ ...
repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
% PREPARE FILTER
[H, Hs] = filter(nelx, nely, nelz, rmin);
% INITIALIZE ITERATION
x = repmat(volfrac, [nely, nelx, nelz]); y = repmat(0.5, [nely, nelx, nelz]);
xPhys = x; yPhys = y;
loop = 0; change = 1;
rstPath = '.\rstData\';
% Initialize MMA Optimizer
m = 1; n = 2*nele;
xmin = [0.1*ones(nele,1); zeros(nele,1)]; % Column vector with the lower bounds for variables x_j.
xmax = [0.6*ones(nele,1); ones(nele,1)]; % Column vector with the upper bounds for variables x_j.
xval = [xPhys(:); yPhys(:)]; % design variables
xold1 = xval(:); xold2 = xold1(:);
low = ones(n,1); upp = ones(n,1);
a0 = 1; a_mma = zeros(m,1); c_mma = 5000*ones(m,1); d_mma = zeros(m,1);
% START ITERATION
while change > tolx && loop < maxloop
loop = loop+1;
% FE-ANALYSIS AND SENSITIVITY ANALYSIS
[K_cell, K_dx1_cell, K_dx2_cell] = arrayfun(@KE_matrix, xPhys(:)', yPhys(:)', 'un', 0);
KALL = reshape(1/8*cell2mat(K_cell), [24*24, nele]);
sK = reshape(KALL, 24*24*nele, 1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
F_nodes = 1/8*cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_cell,'un', 0)');
ce = reshape(sum(F_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
c = sum(sum(sum(ce)));

F_dx1_nodes = 1/8*cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dx1_cell,'un', 0)');
ce_dx1 = reshape(sum(F_dx1_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
dc_x = -ce_dx1;

F_dx2_nodes = 1/8*cell2mat(cellfun(@(x,y)x*y, num2cell(U(edofMat),2)', K_dx2_cell,'un', 0)');
ce_dx2 = reshape(sum(F_dx2_nodes.*U(edofMat), 2), [nely, nelx, nelz]);
dc_y = -ce_dx2;
dv_x = ones(nely,nelx,nelz);
% FILTERING AND MODIFICATION OF SENSITIVITIES
dc_x(:) = H*(x(:).*dc_x(:))./Hs./max(1e-3,x(:));
f0val = c; df0dx = [dc_x(:); dc_y(:)]; df0dx2 = zeros(n, 1);
fval = sum(xPhys(:))/(volfrac*nele) - 1;
dfdx = [dv_x(:)'/(volfrac*nele), zeros(1,nele)]; dfdx2 = 0*dfdx;
[xmma,~,~,~,~,~,~,~,~,low,upp] = ...
mmasub(m,n,loop,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a_mma,c_mma,d_mma);
xold2 = xold1; xold1 = xval; change = max(abs(xmma-xval)); xval = xmma;
xnew = reshape(xval(1:nele), nely, nelx, nelz);
ynew = reshape(xval(nele+1:end), nely, nelx, nelz);
% FILTERING AND MODIFICATION OF VARIABLES
ynew(:) = H*(ynew(:)./Hs); ynew(ynew > 1.0) = 1.0;
x = xnew; y = ynew;
% PRINT RESULTS
fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f \n',loop,c,mean(xPhys(:)),change);
xPhys = x; yPhys = y;
end
% SAVE THE OPTIMAL RESULTS
save('rst.mat', 'x', 'y', 'nelx', 'nely', 'nelz');
end
890 C. Wang et al.% PREPARE FILTER
function [H, Hs] = filter(nelx, nely, nelz, rmin)
nele = nelx*nely*nelz;
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
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
end
% CALCULATION OF GEOMETRIC MATRIX
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
Bs(:,:,i) = [ a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i), 0, 0;
0, -d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i), 0;
0, 0, g*diff(Ns(i),t)-h*diff(Ns(i), t)+l*diff(Ns(i),n);
-d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i),n), a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i),n), 0;
0, g*diff(Ns(i),s)-h*diff(Ns(i),t)+l*diff(Ns(i),n), -d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i),n);
g*diff(Ns(i),s)-h*diff(Ns(i),t)+l*diff(Ns(i),n), 0,
a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i),n)]/Jdet;
end
Concurrent design of hierarchical structures with three-dimensional parameterized lattice microstructures... 891B = [Bs(:,:,1),Bs(:,:,2),Bs(:,:,3),Bs(:,:,4),Bs(:,:,5),Bs(:,:,6),Bs(:,:,7),Bs(:,:,8)];
B = matlabFunction(B);
end
% CALCULATION OF STIFFNESS MATRIX AND DERIVATIVE
% The stiffness matrix depends on the topologies of the lattice microstructures and the PILM model
function [KE, KE_dx1, KE_dx2] = KE_matrix(x1, x2)
global B;
d11 = -0.172 + 17.4*x2 + 0.6326*x1 -129.3*x2^2 -10.26*x2*x1 + 166.2*x1^2 + 344*x2^3 + 146.5*x2^2*x1
+ ...
38.33*x2*x1^2 -426.1*x1^3 -374.2*x2^4 -273.3*x2^3*x1 + 35.33*x2^2*x1^2 - 83.48*x2*x1^3 +
646.9*x1^4 + ...
141.9*x2^5 + 169.3*x2^4*x1 -75.7*x2^3*x1^2 + 42.17*x2^2*x1^3 + 10.92*x2*x1^4 -252.1*x1^5;
d44 = -0.08025 + 3.865*x2 + 8.943*x1 - 23.02*x2^2 - 36.99*x2*x1 + 61.1*x1^2 + 56.08*x2^3 +
84.63*x2^2*x1 +...
18.81*x2*x1^2 - 154.5*x1^3 - 58.9*x2^4 - 108.4*x2^3*x1 - 68.18*x2^2*x1^2 + 70.76*x2*x1^3 +
202.9*x1^4 + 21.84*x2^5 +...
55.7*x2^4*x1 - 4.93*x2^3*x1^2 + 56.04*x2^2*x1^3 - 67.2*x2*x1^4 - 79.76*x1^5;
d12 = -0.132 + 10.59*x2 + 2.776*x1 -64.49*x2^2 -59.67*x2*x1 + 129.8*x1^2 + 152.5*x2^3 +
179.3*x2^2*x1 -...
2.775*x2*x1^2 - 329.8*x1^3 - 154.4*x2^4 - 230.6*x2^3*x1 - 100.2*x2^2*x1^2 + 155.9*x2*x1^3 +
325.1*x1^4 +...
55.57*x2^5 + 115*x2^4*x1 -5.98*x2^3*x1^2 + 72.37*x2^2*x1^3 - 123*x2*x1^4 - 69.87*x1^5;
d11_dx1 = 0.6326 - 10.26*x2 + 2*166.2*x1 + 146.5*x2^2 + 2*38.33*x2*x1 - 3*426.1*x1^2 - 273.3*x2^3 +
2*35.33*x2^2*x1 -...
3*83.48*x2*x1^2 + 4*646.9*x1^3 + 169.3*x2^4 - 2*75.7*x2^3*x1 + 3*42.17*x2^2*x1^2 +
4*10.92*x2*x1^3 -5*252.1*x1^4;
d44_dx1 = 8.943 - 36.99*x2 + 2*61.1*x1 + 84.63*x2^2+ 2*18.81*x2*x1 - 3*154.5*x1^2 - 108.4*x2^3 -
2*68.18*x2^2*x1 +...
3*70.76*x2*x1^2 + 4*202.9*x1^3 + 55.7*x2^4 - 2*4.93*x2^3*x1 + 3*56.04*x2^2*x1^2 -
4*67.2*x2*x1^3 - 5*79.76*x1^4;
d12_dx1 = 2.776 - 59.67*x2 + 2*129.8*x1 + 179.3*x2^2 - 2*2.775*x2*x1 - 3*329.8*x1^2 - 230.6*x2^3 -
2*100.2*x2^2*x1 +...
3*155.9*x2*x1^2 + 4*325.1*x1^3 + 115*x2^4 - 2*5.98*x2^3*x1 + 3*72.37*x2^2*x1^2 -
4*123*x2*x1^3 - 5*69.87*x1^4;
d11_dx2 = 17.4 - 2*129.3*x2 - 10.26*x2 + 3*344*x2^2 + 2*146.5*x2*x1 + 38.33*x1^2 - 4*374.2*x2^3 -
3*273.3*x2^2*x1 +...
2*35.33*x2*x1^2 - 83.48*x1^3 + 5*141.9*x2^4 + 4*169.3*x2^3*x1 - 3*75.7*x2^2*x1^2 +
2*42.17*x2*x1^3 + 10.92*x1^4;
d44_dx2 = 3.865 - 2*23.02*x2 - 36.99*x1 + 3*56.08*x2^2 + 2*84.63*x2*x1 + 18.81*x1^2 - 4*58.9*x2^3 -
3*108.4*x2^2*x1 -...
2*68.18*x2*x1^2 + 70.76*x1^3 + 5*21.84*x2^4 + 4*55.7*x2^3*x1 - 3*4.93*x2^2*x1^2 +
2*56.04*x2*x1^3 - 67.2*x1^4;
d12_dx2 = 10.59 - 2*64.49*x2 - 59.67*x1 + 3*152.5*x2^2 + 2*179.3*x2*x1 - 2.775*x1^2 - 4*154.4*x2^3 -
3*230.6*x2^2*x1 -...
2*100.2*x2*x1^2 + 155.9*x1^3 + 5*55.57*x2^4 + 4*115*x2^3*x1 - 3*5.98*x2^2*x1^2 +
2*72.37*x2*x1^3 - 123*x1^4;
D = [d11, d12, d12, 0, 0, 0;
d12, d11, d12, 0, 0, 0;
d12, d12, d11, 0, 0, 0;
0, 0, 0, d44, 0, 0;
0, 0, 0, 0, d44, 0;
0, 0, 0, 0, 0, d44];
D_dx1 = [d11_dx1, d12_dx1, d12_dx1, 0, 0, 0;
d12_dx1, d11_dx1, d12_dx1, 0, 0, 0;
d12_dx1, d12_dx1, d11_dx1, 0, 0, 0;
0, 0, 0, d44_dx1, 0, 0;
0, 0, 0, 0, d44_dx1, 0;
0, 0, 0, 0, 0, d44_dx1];
D_dx2 = [d11_dx2, d12_dx2, d12_dx2, 0, 0, 0;
d12_dx2, d11_dx2, d12_dx2, 0, 0, 0;
d12_dx2, d12_dx2, d11_dx2, 0, 0, 0;
0, 0, 0, d44_dx2, 0, 0;
0, 0, 0, 0, d44_dx2, 0;
0, 0, 0, 0, 0, d44_dx2];
x = [-sqrt(15)/5, 0, sqrt(15)/5];
w = [5/9, 8/9, 5/9];
ke = zeros(24);
ke_dx1 = zeros(24);
ke_dx2 = zeros(24);
[i, j, k] = meshgrid(1:3, 1:3, 1:3);
for m = 1:27
ke = ke + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D * B(x(i(m)),
x(j(m)), x(k(m)));
ke_dx1 = ke_dx1 + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D_dx1 * B(x(i(m)), x(j(m)),
x(k(m)));
ke_dx2 = ke_dx2 + w(i(m))*w(j(m))*w(k(m))*B(x(i(m)), x(j(m)), x(k(m)))' * D_dx2 * B(x(i(m)), x(j(m)),
x(k(m)));
end
KE = ke; KE_dx1 = ke_dx1; KE_dx2 = 10*ke_dx2;
end