LHid=fopen('gridNoNoise.txt');
LHpoints = fscanf(LHid,' %f %f %f %f %f %f',[6 inf]);
fclose('all')
initDes=1:24;
locfilename=['objpoint6'];
locfileID=fopen(locfilename,'wt');

for j=501:600
    for i=1:size(initDes,2)
        [~,obji,~]=unitCell8tz_3D(9,9,9,LHpoints(1,j),3,1.5,2,LHpoints(2,j),LHpoints(3,j),LHpoints(4,j),LHpoints(5,j),LHpoints(6,j),initDes(i),0.5); %changer après
        fprintf(locfileID,'%14.10f ',obji);
        i*j
    end
    fprintf(locfileID,'\n');
end

fclose('all');


