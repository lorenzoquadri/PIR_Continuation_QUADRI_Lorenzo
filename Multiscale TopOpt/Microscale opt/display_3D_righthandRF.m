function display_3D_righthandRF(rho)
[nelz,nely,nelx] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
% set(gca, 'XDir','reverse');
% set(gca, 'YDir','reverse');
axis equal; axis tight; axis on; xlabel('x'); ylabel('y'); zlabel('z');
box on; view([135,45]);  % to show axis

for i = 1:nelx
    x = (i-1)*hx;
    for j = 1:nely 
    %for i = round(nely/2):nely % to show half cube
        y = (j-1)*hy;
        for k = 1:nelz
        % for j = round(nelz/2):nelz % to show half cube cut at y/2
            %y = nelz*hy - (j-1)*hy;
            z = (k-1)*hz;
            if (rho(nelz+1-k,j,nelx+1-i) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y+hx z; x+hx y+hx z; x+hx y z; x y z+hx;x y+hx z+hx; x+hx y+hx z+hx;x+hx y z+hx];
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(nelz+1-k,j,nelx+1-i)),0.2+0.8*(1-rho(nelz+1-k,j,nelx+1-i)),0.2+0.8*(1-rho(nelz+1-k,j,nelx+1-i))]);
                hold on;
            end
        end
    end
end

 
end