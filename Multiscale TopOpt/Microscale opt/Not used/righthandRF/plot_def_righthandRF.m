function plot_def_righthandRF(U,nelx,nely,nelz)
if nargin==4
    %
%     nelx=0;
    % 
    
    [y_in,z_in,x_in]=meshgrid(0:nely,nelz:-1:0,0:nelx);

    
    %
%     x_fin(:)=x_in(:)+U(1:3:3*(nely+1)*(nelz+1));
%     y_fin(:)=y_in(:)-U(2:3:3*(nely+1)*(nelz+1));
%     z_fin(:)=z_in(:);
    %
    
    x_fin(:)=x_in(:)+U(1:3:end);
    y_fin(:)=y_in(:)+U(2:3:end);
    z_fin(:)=z_in(:)+U(3:3:end);
    plot3(x_fin(:),y_fin(:),z_fin(:),'ro');
    hold on
    plot3(x_fin(1),y_fin(1),z_fin(1),'r*',x_fin(nelz+1),y_fin(nelz+1),z_fin(nelz+1),'r*',x_fin((nelz+1)*nely+1),y_fin((nelz+1)*nely+1),z_fin((nelz+1)*nely+1),'r*',x_fin((nelz+1)*(nely+1)),y_fin((nelz+1)*(nely+1)),z_fin((nelz+1)*(nely+1)),'r*',...
        x_fin((nelz+1)*(nely+1)*nelx+1),y_fin((nelz+1)*(nely+1)*nelx+1),z_fin((nelz+1)*(nely+1)*nelx+1),'r*',x_fin((nelz+1)*(nely+1)*nelx+nelz+1),y_fin((nelz+1)*(nely+1)*nelx+nelz+1),z_fin((nelz+1)*(nely+1)*nelx+nelz+1),'r*',x_fin((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),y_fin((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),z_fin((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),'r*',x_fin((nelz+1)*(nely+1)*(nelx+1)),y_fin((nelz+1)*(nely+1)*(nelx+1)),z_fin((nelz+1)*(nely+1)*(nelx+1)),'r*');
    
    plot3(x_in(:),y_in(:),z_in(:),'bo');
    plot3(x_in(1),y_in(1),z_in(1),'b*',x_in(nelz+1),y_in(nelz+1),z_in(nelz+1),'b*',x_in((nelz+1)*nely+1),y_in((nelz+1)*nely+1),z_in((nelz+1)*nely+1),'b*',x_in((nelz+1)*(nely+1)),y_in((nelz+1)*(nely+1)),z_in((nelz+1)*(nely+1)),'b*',...
       x_in((nelz+1)*(nely+1)*nelx+1),y_in((nelz+1)*(nely+1)*nelx+1),z_in((nelz+1)*(nely+1)*nelx+1),'b*',x_in((nelz+1)*(nely+1)*nelx+nelz+1),y_in((nelz+1)*(nely+1)*nelx+nelz+1),z_in((nelz+1)*(nely+1)*nelx+nelz+1),'b*',x_in((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),y_in((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),z_in((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),'b*',x_in((nelz+1)*(nely+1)*(nelx+1)),y_in((nelz+1)*(nely+1)*(nelx+1)),z_in((nelz+1)*(nely+1)*(nelx+1)),'b*');
    
    set(gcf,'Name','ISO display','NumberTitle','off');
    axis equal; axis tight; axis on; xlabel('x'); ylabel('y'); zlabel('z');
    box on; view([120,30]);  % to show axis
elseif nargin == 3
    nelx=0;
    [x_in,y_in,z_in]=meshgrid(0:nelz,0:nely,nelx);

    plot3(x_in(:),z_in(:),y_in(:),'bo');
    hold on
    plot3(x_in(1),z_in(1),y_in(1),'b*',x_in(nelz+1),z_in(nelz+1),y_in(nelz+1),'b*',x_in((nelz+1)*nely+1),z_in((nelz+1)*nely+1),y_in((nelz+1)*nely+1),'b*',x_in((nelz+1)*(nely+1)),z_in((nelz+1)*(nely+1)),y_in((nelz+1)*(nely+1)),'b*',...
        x_in((nelz+1)*(nely+1)*nelx+1),z_in((nelz+1)*(nely+1)*nelx+1),y_in((nelz+1)*(nely+1)*nelx+1),'b*',x_in((nelz+1)*(nely+1)*nelx+nelz+1),z_in((nelz+1)*(nely+1)*nelx+nelz+1),y_in((nelz+1)*(nely+1)*nelx+nelz+1),'b*',x_in((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),z_in((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),y_in((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),'b*',x_in((nelz+1)*(nely+1)*(nelx+1)),z_in((nelz+1)*(nely+1)*(nelx+1)),y_in((nelz+1)*(nely+1)*(nelx+1)),'b*');

    x_fin(:)=x_in(:)+U(1:2:end);
    y_fin(:)=y_in(:)-U(2:2:end);
    z_fin(:)=z_in(:);
    plot3(x_fin(:),z_fin(:),y_fin(:),'ro');
    plot3(x_fin(1),z_fin(1),y_fin(1),'r*',x_fin(nelz+1),z_fin(nelz+1),y_fin(nelz+1),'r*',x_fin((nelz+1)*nely+1),z_fin((nelz+1)*nely+1),y_fin((nelz+1)*nely+1),'r*',x_fin((nelz+1)*(nely+1)),z_fin((nelz+1)*(nely+1)),y_fin((nelz+1)*(nely+1)),'r*',...
        x_fin((nelz+1)*(nely+1)*nelx+1),z_fin((nelz+1)*(nely+1)*nelx+1),y_fin((nelz+1)*(nely+1)*nelx+1),'r*',x_fin((nelz+1)*(nely+1)*nelx+nelz+1),z_fin((nelz+1)*(nely+1)*nelx+nelz+1),y_fin((nelz+1)*(nely+1)*nelx+nelz+1),'r*',x_fin((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),z_fin((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),y_fin((nelz+1)*(nely+1)*nelx+(nelz+1)*nely+1),'r*',x_fin((nelz+1)*(nely+1)*(nelx+1)),z_fin((nelz+1)*(nely+1)*(nelx+1)),y_fin((nelz+1)*(nely+1)*(nelx+1)),'r*');

    set(gcf,'Name','ISO display','NumberTitle','off');
    axis equal; axis tight; axis on; xlabel('x'); ylabel('y'); zlabel('z');
    box on; view([120,30]);  % to show axis
end