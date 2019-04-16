% AMRLevel = 4;
numCPU = 24;
numTimeSteps = 1001;

fileID = fopen('collectionB64.pvd','w');

fprintf(fileID, '<?xml version="1.0"?>\n');
fprintf(fileID, '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n');
fprintf(fileID,'    <Collection>\n');
for i = 1:numTimeSteps
%     for j = 1:AMRLevel+1
%         fprintf(fileID,num2str(j-1));
%         fprintf(fileID,'_TIME');
%         fprintf(fileID,num2str(i-1));
%         fprintf(fileID,'.vtu');
%         fprintf(fileID,'\n');
    if (mod((i-1),2) == 0)
        for j = 1:numCPU
            fprintf(fileID,'        <DataSet timestep="');
            fprintf(fileID,int2str(i-1));
            fprintf(fileID,'" part="');
            fprintf(fileID,int2str(j-1));
            fprintf(fileID,'" file="resultB64_CPU');
            fprintf(fileID,int2str(j-1));
            fprintf(fileID,'_AMR0_TIME');
            fprintf(fileID,int2str(i-1));
            fprintf(fileID,'.vtr"/>');
            fprintf(fileID,'\n');
        end
    end
end
fprintf(fileID,'    </Collection>\n');
fprintf(fileID,'</VTKFile>\n');
fclose(fileID);