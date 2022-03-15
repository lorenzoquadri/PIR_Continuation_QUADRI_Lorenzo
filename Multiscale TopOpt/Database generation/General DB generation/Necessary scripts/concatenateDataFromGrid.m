locfilename='databaseGridNoNoise.txt';
locfilename2='objGridNoNoise.txt';
locfilename3='struGridNoNoise.txt';

for i=1:61200
    i
    fileID=fopen(locfilename,'a+');
    fileID2=fopen(locfilename2,'a+');
    fileID3=fopen(locfilename3,'a+');
    
    pointfilename=['Results/point',num2str(i)];
    pointfileID=fopen(pointfilename);
    pointdata = fscanf(pointfileID,'%f',[27 inf]);
    
    pointobjfilename=['Results/objpoint',num2str(i)];
    pointobjfileID=fopen(pointobjfilename);
    pointobjdata = fscanf(pointobjfileID,'%f',[12 inf]);
    
    pointstrufilename=['Results/strupoint',num2str(i)];
    pointstrufileID=fopen(pointstrufilename);
    pointstrudata = fscanf(pointstrufileID,'%f',[735 inf]);
    
    fprintf(fileID,'%14.10f',pointdata');
    fprintf(fileID,'\n');
    fprintf(fileID2,'%14.10f',pointobjdata');
    fprintf(fileID2,'\n');
    fprintf(fileID3,'%14.10f',pointstrudata');
    fprintf(fileID3,'\n');
    
    fclose('all');
end




