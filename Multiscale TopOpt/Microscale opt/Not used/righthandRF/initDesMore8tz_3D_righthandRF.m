function x=initDesMore8tz_3D_righthandRF(nelx,nely,nelz,volfrac,initDes,thickInit)
if initDes==1
    x=ones(nelz,nely,nelx);
elseif initDes==2
    x = volfrac*ones(nelz,nely,nelx);
elseif initDes==3
    x = 0.5*ones(nelz,nely,nelx);
% all edges and zeros
elseif initDes==4
    x=zeros(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
            end
        end
    end
% all edges and volfrac
elseif initDes==5
    x = volfrac*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
            end
        end
    end
% all edges and 0.5
elseif initDes==6
    x = 0.5*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
            end
        end
    end
% all 2D diagonals and zeros
elseif initDes==7
    x=zeros(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
            end
        end
    end   
% all 2D diagonals and volfrac
elseif initDes==8
    x = volfrac*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
            end
        end
    end   
% all 2D diagonals and 0.5
elseif initDes==9
    x = 0.5*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
            end
        end
    end 
% all 3D diagonals and 0
elseif initDes==10
    x=zeros(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end  
% all 3D diagonals and volfrac
elseif initDes==11
    x = volfrac*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end
% all 3D diagonals and 0.5
elseif initDes==12
    x = 0.5*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end
% all edges and 2D diagonals and zeros
elseif initDes==13
    x=zeros(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
                % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
            end
        end
    end
% all edges and 2D diagonals and volfrac
elseif initDes==14
    x = volfrac*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
                % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
            end
        end
    end
% all edges and 2D diagonals and 0.5
elseif initDes==15
    x = 0.5*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
                % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
            end
        end
    end
% all edges and 3D diagonals and zeros
elseif initDes==16
    x=zeros(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end
% all edges and and 3D diagonals volfrac
elseif initDes==17
    x = volfrac*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end
% all edges and and 3D diagonals 0.5
elseif initDes==18
    x = 0.5*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
                % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end
% all 2D diagonals and 3D diagonals and zeros
elseif initDes==19
    x=zeros(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
               
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end   
% all 2D diagonals  and 3D diagonals and volfrac
elseif initDes==20
    x = volfrac*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
               
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end   
% all 2D diagonals  and 3D diagonals and 0.5
elseif initDes==21
    x = 0.5*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
               
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
            end
        end
    end 
% all 2D diagonals and 3D diagonals and edges zeros
elseif initDes==22
    x=zeros(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
               
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
            end
        end
    end   
% all 2D diagonals  and 3D diagonals and edges volfrac
elseif initDes==23
    x = volfrac*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
               
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
            end
        end
    end   
% all 2D diagonals  and 3D diagonals and edges and 0.5
elseif initDes==24
    x = 0.5*ones(nelz,nely,nelx);
    for ely=1:nely
        for elz=1:nelz
            for elx=1:nelx
               % xy sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xy sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides downwards diagonals
               if ((ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && elz<=thickInit) || (ely-elx < thickInit/sqrt(2) && ely-elx > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % xz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && elz<=thickInit) || ((ely-1)/(nely-1)+(elx-1)/(nelx-1) < 1+ 2/(nely+nelx)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elx-1)/(nelx-1) > 1- 2/(nely+nelx)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(elz,ely,elx)=1;
               end
               % yz sides downwards diagonals
               if ((elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && ely<=thickInit) || (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2) && nelz-ely<thickInit))
                x(elz,ely,elx)=1;
               end
                % yz sides upwards diagonals
               if (((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && ely<=thickInit) || ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(elz,ely,elx)=1;
               end
               
                % xy downwards + yz downwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz downwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && (elz-elx < thickInit/sqrt(2) && elz-elx > - thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               
               % xy downwards + yz upwards
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2))))
                x(elz,ely,elx)=1;
               end
               % xy upwards + yz upwards
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)) && ((elz-1)/(nelz-1)+(elx-1)/(nelx-1) < 1+ 2/(nelz+nelx)*thickInit/sqrt(2) && (elz-1)/(nelz-1)+(elx-1)/(nelx-1) > 1- 2/(nelz+nelx)*thickInit/sqrt(2)))
                x(elz,ely,elx)=1;
               end
               % y-axis edges
                if((ely<=thickInit && elx<=thickInit) || (nely-ely<thickInit && nelx-elx<thickInit) || (nely-ely<thickInit && elx<=thickInit) || (ely<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % x-axis edges
                if((elz<=thickInit && elx<=thickInit) || (nelz-elz<thickInit && nelx-elx<thickInit) || (nelz-elz<thickInit && elx<=thickInit) || (elz<=thickInit && nelx-elx<thickInit))
                    x(elz,ely,elx)=1;
                end
                % z-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(elz,ely,elx)=1;
                end
            end
        end
    end 
% test case
elseif initDes==25
    x = zeros(nelz,nely,nelx);
    x(1,:,1)=1;
    x(nelz,:,1)=1;
    x(nelz,:,nelx)=1;
    x(1,:,nelx)=1;
    % test case: to prove Poisson effect
elseif initDes==26
    x=ones(nelz,nely,nelx);
    for i=1:nelz
        for j=1:nely
            for k=2:nelx-1
                if (i==1||j==1||i==nelz||j==nely)
                    x(i,j,k)=0;
                end
            end
        end
    end
else
    error('initDes should be in 1..24')
end
figure(1)
display_3D_righthandRF(x)
