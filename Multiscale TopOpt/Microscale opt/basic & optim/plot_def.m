function plot_def(U,nelx,nely,nelz)
if nargin==4
    %
%     nelz=0;
    % 
    
    
    [x_in,y_in,z_in]=meshgrid(0:nelx,nely:-1:0,0:nelz);
    
    %
%     x_fin(:)=x_in(:)+U(1:3:3*(nelx+1)*(nely+1));
%     y_fin(:)=y_in(:)-U(2:3:3*(nelx+1)*(nely+1));
%     z_fin(:)=z_in(:);
    %
    
    x_fin(:)=x_in(:)+U(1:3:end);
    y_fin(:)=y_in(:)+U(2:3:end);
    z_fin(:)=z_in(:)+U(3:3:end);
    plot3(x_fin(:),z_fin(:),y_fin(:),'ro');
    hold on
    plot3(x_fin(1),z_fin(1),y_fin(1),'r*',x_fin(nely+1),z_fin(nely+1),y_fin(nely+1),'r*',x_fin((nely+1)*nelx+1),z_fin((nely+1)*nelx+1),y_fin((nely+1)*nelx+1),'r*',x_fin((nely+1)*(nelx+1)),z_fin((nely+1)*(nelx+1)),y_fin((nely+1)*(nelx+1)),'r*',...
        x_fin((nely+1)*(nelx+1)*nelz+1),z_fin((nely+1)*(nelx+1)*nelz+1),y_fin((nely+1)*(nelx+1)*nelz+1),'r*',x_fin((nely+1)*(nelx+1)*nelz+nely+1),z_fin((nely+1)*(nelx+1)*nelz+nely+1),y_fin((nely+1)*(nelx+1)*nelz+nely+1),'r*',x_fin((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),z_fin((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),y_fin((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),'r*',x_fin((nely+1)*(nelx+1)*(nelz+1)),z_fin((nely+1)*(nelx+1)*(nelz+1)),y_fin((nely+1)*(nelx+1)*(nelz+1)),'r*');
    
    plot3(x_in(:),z_in(:),y_in(:),'bo');
    plot3(x_in(1),z_in(1),y_in(1),'b*',x_in(nely+1),z_in(nely+1),y_in(nely+1),'b*',x_in((nely+1)*nelx+1),z_in((nely+1)*nelx+1),y_in((nely+1)*nelx+1),'b*',x_in((nely+1)*(nelx+1)),z_in((nely+1)*(nelx+1)),y_in((nely+1)*(nelx+1)),'b*',...
        x_in((nely+1)*(nelx+1)*nelz+1),z_in((nely+1)*(nelx+1)*nelz+1),y_in((nely+1)*(nelx+1)*nelz+1),'b*',x_in((nely+1)*(nelx+1)*nelz+nely+1),z_in((nely+1)*(nelx+1)*nelz+nely+1),y_in((nely+1)*(nelx+1)*nelz+nely+1),'b*',x_in((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),z_in((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),y_in((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),'b*',x_in((nely+1)*(nelx+1)*(nelz+1)),z_in((nely+1)*(nelx+1)*(nelz+1)),y_in((nely+1)*(nelx+1)*(nelz+1)),'b*');
    
    set(gcf,'Name','ISO display','NumberTitle','off');
    % set(gca, 'ZDir','reverse');

    axis equal; axis tight; axis on; xlabel('x'); ylabel('z'); zlabel('y');
    box on; view([30,30]); pause(1e-6); % to show axis % view([0,0]);
elseif nargin == 3
    nelz=0;
    [x_in,y_in,z_in]=meshgrid(0:nely,0:nelx,nelz);

    plot3(x_in(:),z_in(:),y_in(:),'bo');
    hold on
    plot3(x_in(1),z_in(1),y_in(1),'b*',x_in(nely+1),z_in(nely+1),y_in(nely+1),'b*',x_in((nely+1)*nelx+1),z_in((nely+1)*nelx+1),y_in((nely+1)*nelx+1),'b*',x_in((nely+1)*(nelx+1)),z_in((nely+1)*(nelx+1)),y_in((nely+1)*(nelx+1)),'b*',...
        x_in((nely+1)*(nelx+1)*nelz+1),z_in((nely+1)*(nelx+1)*nelz+1),y_in((nely+1)*(nelx+1)*nelz+1),'b*',x_in((nely+1)*(nelx+1)*nelz+nely+1),z_in((nely+1)*(nelx+1)*nelz+nely+1),y_in((nely+1)*(nelx+1)*nelz+nely+1),'b*',x_in((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),z_in((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),y_in((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),'b*',x_in((nely+1)*(nelx+1)*(nelz+1)),z_in((nely+1)*(nelx+1)*(nelz+1)),y_in((nely+1)*(nelx+1)*(nelz+1)),'b*');

    x_fin(:)=x_in(:)+U(1:2:end);
    y_fin(:)=y_in(:)-U(2:2:end);
    z_fin(:)=z_in(:);
    plot3(x_fin(:),z_fin(:),y_fin(:),'ro');
    plot3(x_fin(1),z_fin(1),y_fin(1),'r*',x_fin(nely+1),z_fin(nely+1),y_fin(nely+1),'r*',x_fin((nely+1)*nelx+1),z_fin((nely+1)*nelx+1),y_fin((nely+1)*nelx+1),'r*',x_fin((nely+1)*(nelx+1)),z_fin((nely+1)*(nelx+1)),y_fin((nely+1)*(nelx+1)),'r*',...
        x_fin((nely+1)*(nelx+1)*nelz+1),z_fin((nely+1)*(nelx+1)*nelz+1),y_fin((nely+1)*(nelx+1)*nelz+1),'r*',x_fin((nely+1)*(nelx+1)*nelz+nely+1),z_fin((nely+1)*(nelx+1)*nelz+nely+1),y_fin((nely+1)*(nelx+1)*nelz+nely+1),'r*',x_fin((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),z_fin((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),y_fin((nely+1)*(nelx+1)*nelz+(nely+1)*nelx+1),'r*',x_fin((nely+1)*(nelx+1)*(nelz+1)),z_fin((nely+1)*(nelx+1)*(nelz+1)),y_fin((nely+1)*(nelx+1)*(nelz+1)),'r*');

    set(gcf,'Name','ISO display','NumberTitle','off');
    set(gca, 'ZDir','reverse');

    axis equal; axis tight; axis on; xlabel('x'); ylabel('z'); zlabel('y');
    box on; view([30,30]); pause(1e-6); % to show axis
end