function verify_struct
for i=2500:4096
    i
    filename=['strupoint',num2str(i)];
    fileID=fopen(filename);
    data = fscanf(fileID,'%f',[1 inf]);
    microv=data(7:end);
    micro=reshape(microv,9,9,9);
    display_3D(micro)
    
end
end