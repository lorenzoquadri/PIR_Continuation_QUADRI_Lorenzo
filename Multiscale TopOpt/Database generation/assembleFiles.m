nparts=16;

objtfileID=fopen('objGridNoNoisetot.txt','w');
for i=1:nparts
    objfilename=['DB2Assemble/objGridNoNoise',num2str(i),'.txt'];
    objfileID=fopen(objfilename);
    objdata = fscanf(objfileID,'%f',[30 inf]);
    objdata=objdata';
    size(objdata,1)
    for j=1:size(objdata,1)
        fprintf(objtfileID,'%14.10f',objdata(j,:))
        fprintf(objtfileID,'\n')
    end
end

dbtfileID=fopen('databaseGridNoNoisetot.txt','w');
for i=1:nparts
    dbfilename=['DB2Assemble/databaseGridNoNoise',num2str(i),'.txt'];
    dbfileID=fopen(dbfilename);
    dbdata = fscanf(dbfileID,'%f',[27 inf]);
    dbdata=dbdata';
    for j=1:size(dbdata,1)
        fprintf(dbtfileID,'%14.10f',dbdata(j,:))
        fprintf(dbtfileID,'\n')
    end
end

strutfileID=fopen('struGridNoNoisetot.txt','w');
for i=1:nparts
    strufilename=['DB2Assemble/struGridNoNoise',num2str(i),'.txt'];
    strufileID=fopen(strufilename);
    strudata = fscanf(strufileID,'%f',[735 inf]);
    strudata=strudata';
    for j=1:size(strudata,1)
        fprintf(strutfileID,'%14.10f',strudata(j,:))
        fprintf(strutfileID,'\n')
    end
end