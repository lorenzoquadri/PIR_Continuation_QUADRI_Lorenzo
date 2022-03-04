LHid=fopen('gridNoNoise.txt');
LHpoints = fscanf(LHid,'%f %f %f %f %f %f',[6 inf]);
fclose('all');
initDes=[1,6,9,14,18,21]; % Select best initDes
parfor j=14001:16000
% for j=1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j
    locfilename=['point',num2str(j)];
    locfilename2=['objpoint',num2str(j)];
    locfilename3=['strupoint',num2str(j)];
    locfileID=fopen(locfilename,'w');
    locfileID2=fopen(locfilename2,'w');
    locfileID3=fopen(locfilename3,'w');
    obj=0;
    Qopt=zeros(1,21);
    Sopt=zeros(1,729); %changer après
    fprintf(locfileID2,'%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f ',LHpoints(1,j),LHpoints(2,j),LHpoints(3,j),LHpoints(4,j),LHpoints(5,j),LHpoints(6,j));
    for i=1:size(initDes,2)
        [Q,obji,microstruct]=unitCell8tz_3D(9,9,9,LHpoints(1,j),3,1.5,2,LHpoints(2,j),LHpoints(3,j),LHpoints(4,j),LHpoints(5,j),LHpoints(6,j),initDes(i),0.5); %changer après
        fprintf(locfileID2,'%14.10f ',obji);
        if obji<obj
            obj=obji;
            Qopt=Q;
            Sopt=microstruct(:)';
        end
    end
    fprintf(locfileID,'%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f',...
    LHpoints(1,j),LHpoints(2,j),LHpoints(3,j),LHpoints(4,j),LHpoints(5,j),LHpoints(6,j),...
    Qopt(1,1),Qopt(1,2),Qopt(1,3),Qopt(1,4),Qopt(1,5),Qopt(1,6),Qopt(2,2),Qopt(2,3),Qopt(2,4),Qopt(2,5),Qopt(2,6),Qopt(3,3),Qopt(3,4),Qopt(3,5),Qopt(3,6),Qopt(4,4),Qopt(4,5),Qopt(4,6),Qopt(5,5),Qopt(5,6),Qopt(6,6));
    fprintf(locfileID3,'%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f',LHpoints(1,j),LHpoints(2,j),LHpoints(3,j),LHpoints(4,j),LHpoints(5,j),LHpoints(6,j));
    fprintf(locfileID3,'%14.10f ',Sopt);
end

% fileID=fopen('databaseGridNoNoise.txt','w');
% fileID2=fopen('objGridNoNoise.txt','w');
% fileID3=fopen('struGridNoNoise.txt','w');
% for i=1:6554%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %for i=1:size(LHpoints,2)
% %for i=14500:14500
%     pointfilename=['point',num2str(i)];
%     pointfileID=fopen(pointfilename);
%     pointdata = fscanf(pointfileID,'%f %f %f %f %f %f %f %f %f',[9 inf]);
%     pointobjfilename=['objpoint',num2str(i)];
%     pointobjfileID=fopen(pointobjfilename);
%     pointobjdata = fscanf(pointobjfileID,'%f %f %f %f %f %f %f %f %f',[11 inf]);
%     pointstrufilename=['strupoint',num2str(i)];
%     pointstrufileID=fopen(pointstrufilename);
%     pointstrudata = fscanf(pointstrufileID,'%f',[10003 inf]);
%     fprintf(fileID,'%14.10f',pointdata')
%     fprintf(fileID,'\n')
%     fprintf(fileID2,'%14.10f',pointobjdata')
%     fprintf(fileID2,'\n')
%     fprintf(fileID3,'%14.10f',pointstrudata')
%     fprintf(fileID3,'\n')
% end


fclose('all');


