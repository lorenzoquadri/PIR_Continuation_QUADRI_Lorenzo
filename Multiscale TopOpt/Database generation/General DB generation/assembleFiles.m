nparts=5

objtfileID=fopen('objGridNoNoisetot.txt','w');
for i=1:nparts
    objfilename=['objGridNoNoise',num2str(i),'.txt'];
    objfileID=fopen(objfilename);
    objdata = fscanf(objfileID,'%f',[45 inf]);
    objdata=objdata';
    for j=1:size(objdata,1)
        fprintf(objtfileID,'%14.10f',objdata(j,:))
        fprintf(objtfileID,'\n')
    end
end

dbtfileID=fopen('databaseGridNoNoisetot.txt','w');
for i=1:nparts
    dbfilename=['databaseGridNoNoise',num2str(i),'.txt'];
    dbfileID=fopen(dbfilename);
    dbdata = fscanf(dbfileID,'%f',[9 inf]);
    dbdata=dbdata';
    for j=1:size(dbdata,1)
        fprintf(dbtfileID,'%14.10f',dbdata(j,:))
        fprintf(dbtfileID,'\n')
    end
end

strutfileID=fopen('struGridNoNoisetot.txt','w');
for i=1:nparts
    strufilename=['struGridNoNoise',num2str(i),'.txt'];
    strufileID=fopen(strufilename);
    strudata = fscanf(strufileID,'%f',[9 inf]);
    strudata=strudata';
    for j=1:size(strudata,1)
        fprintf(strutfileID,'%14.10f',strudata(j,:))
        fprintf(strutfileID,'\n')
    end
end