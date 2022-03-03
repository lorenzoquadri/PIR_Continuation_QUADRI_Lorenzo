%x=initDesMore8tz_3D(nelx,nely,volfrac,initDes,(nelx+nely)/101);
function x=initDesMore8tz_3D(nelx,nely,nelz,volfrac,initDes,thickInit)
if initDes==1
    x=ones(nely,nelx,nelz);
elseif initDes==2
    x = volfrac*ones(nely,nelx,nelz);
elseif initDes==3
    x = 0.5*ones(nely,nelx,nelz);
% all edges and zeros
elseif initDes==4
    x=zeros(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
            end
        end
    end
% all edges and volfrac
elseif initDes==5
    x = volfrac*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
            end
        end
    end
% all edges and 0.5
elseif initDes==6
    x = 0.5*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
            end
        end
    end
% all 2D diagonals and zeros
elseif initDes==7
    x=zeros(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
            end
        end
    end   
% all 2D diagonals and volfrac
elseif initDes==8
    x = volfrac*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
            end
        end
    end   
% all 2D diagonals and 0.5
elseif initDes==9
    x = 0.5*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
            end
        end
    end 
% all 3D diagonals and 0
elseif initDes==10
    x=zeros(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end  
% all 3D diagonals and volfrac
elseif initDes==11
    x = volfrac*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end
% all 3D diagonals and 0.5
elseif initDes==12
    x = 0.5*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end
% all edges and 2D diagonals and zeros
elseif initDes==13
    x=zeros(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
                % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
            end
        end
    end
% all edges and 2D diagonals and volfrac
elseif initDes==14
    x = volfrac*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
                % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
            end
        end
    end
% all edges and 2D diagonals and 0.5
elseif initDes==15
    x = 0.5*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
                % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
            end
        end
    end
% all edges and 3D diagonals and zeros
elseif initDes==16
    x=zeros(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end
% all edges and and 3D diagonals volfrac
elseif initDes==17
    x = volfrac*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end
% all edges and and 3D diagonals 0.5
elseif initDes==18
    x = 0.5*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
                % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end
% all 2D diagonals and 3D diagonals and zeros
elseif initDes==19
    x=zeros(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
               
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end   
% all 2D diagonals  and 3D diagonals and volfrac
elseif initDes==20
    x = volfrac*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
               
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end   
% all 2D diagonals  and 3D diagonals and 0.5
elseif initDes==21
    x = 0.5*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
               
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
            end
        end
    end 
% all 2D diagonals and 3D diagonals and edges zeros
elseif initDes==22
    x=zeros(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
               
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
            end
        end
    end   
% all 2D diagonals  and 3D diagonals and edges volfrac
elseif initDes==23
    x = volfrac*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
               
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
            end
        end
    end   
% all 2D diagonals  and 3D diagonals and edges and 0.5
elseif initDes==24
    x = 0.5*ones(nely,nelx,nelz);
    for elx=1:nelx
        for ely=1:nely
            for elz=1:nelz
               % xy sides downwards diagonals
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && elz<=thickInit) || (elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xy sides upwards diagonals
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && elz<=thickInit) || ((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2) && nelz-elz<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides downwards diagonals
               if ((elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && ely<=thickInit) || (elx-elz < thickInit/sqrt(2) && elx-elz > - thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % xz sides upwards diagonals
               if (((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && ely<=thickInit) || ((elx-1)/(nelx-1)+(elz-1)/(nelz-1) < 1+ 2/(nelx+nelz)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(elz-1)/(nelz-1) > 1- 2/(nelx+nelz)*thickInit/sqrt(2) && nely-ely<thickInit))
                x(ely,elx,elz)=1;
               end
               % yz sides downwards diagonals
               if ((ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && elx<=thickInit) || (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2) && nely-elx<thickInit))
                x(ely,elx,elz)=1;
               end
                % yz sides upwards diagonals
               if (((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && elx<=thickInit) || ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2) && nelx-elx<thickInit))
                x(ely,elx,elz)=1;
               end
               
                % xy downwards + yz downwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz downwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && (ely-elz < thickInit/sqrt(2) && ely-elz > - thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               
               % xy downwards + yz upwards
               if ((elx-ely < thickInit/sqrt(2) && elx-ely > - thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2))))
                x(ely,elx,elz)=1;
               end
               % xy upwards + yz upwards
               if (((elx-1)/(nelx-1)+(ely-1)/(nely-1) < 1+ 2/(nelx+nely)*thickInit/sqrt(2) && (elx-1)/(nelx-1)+(ely-1)/(nely-1) > 1- 2/(nelx+nely)*thickInit/sqrt(2)) && ((ely-1)/(nely-1)+(elz-1)/(nelz-1) < 1+ 2/(nely+nelz)*thickInit/sqrt(2) && (ely-1)/(nely-1)+(elz-1)/(nelz-1) > 1- 2/(nely+nelz)*thickInit/sqrt(2)))
                x(ely,elx,elz)=1;
               end
               % y-axis edges
                if((elx<=thickInit && elz<=thickInit) || (nelx-elx<thickInit && nelz-elz<thickInit) || (nelx-elx<thickInit && elz<=thickInit) || (elx<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % x-axis edges
                if((ely<=thickInit && elz<=thickInit) || (nely-ely<thickInit && nelz-elz<thickInit) || (nely-ely<thickInit && elz<=thickInit) || (ely<=thickInit && nelz-elz<thickInit))
                    x(ely,elx,elz)=1;
                end
                % z-axis edges
                if((elx<=thickInit && ely<=thickInit) || (nelx-elx<thickInit && nely-ely<thickInit) || (nelx-elx<thickInit && ely<=thickInit) || (elx<=thickInit && nely-ely<thickInit))
                    x(ely,elx,elz)=1;
                end
            end
        end
    end 
% test case
elseif initDes==25
    x = zeros(nely,nelx,nelz);
    x(1,:,1)=1;
    x(nely,:,1)=1;
    x(nely,:,nelz)=1;
    x(1,:,nelz)=1;
else
    error('initDes should be in 1..24')
end
%figure(1)
%display_3D(x)
