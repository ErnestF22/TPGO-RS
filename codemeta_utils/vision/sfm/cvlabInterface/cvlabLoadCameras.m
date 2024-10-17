%function [G,P,dim]=cvlabLoadCameras(dirName)
%Load camera pose and projection matrices from files with the "CVLab
%format" in the directory fullfile(dirName,'cameras')
%The poses in G are in the 'reference' interpretation

function [G,P,dim]=cvlabLoadCameras(dirName)
flagNegativeTranslation=false;

listFiles=cvlabGetListFiles(dirName);

NCameras=length(listFiles);

P=zeros(3,4,NCameras);
G=zeros(4,4,NCameras);
G(4,4,:)=1;
dim=zeros(2,NCameras);

for iCamera=1:NCameras
    fileName=fullfile(dirName,'cameras',listFiles{iCamera});
    fid=fopen(fileName,'r');
    if fid<0
        error('Error opening file %s', fileName)
    end
    P(:,:,iCamera)=reshape(fscanf(fid, '%f', 12),3,4);
    G(1:3,1:4,iCamera)=reshape(fscanf(fid, '%f', 12),3,4);
    G(1:3,1:3,iCamera)=G(1:3,1:3,iCamera)';
    if flagNegativeTranslation
        G(1:3,4,iCamera)=-G(1:3,4,iCamera);
    end
    %make sure rotation is orthonormal
    [U,S,V]=svd(G(1:3,1:3,iCamera));
    G(1:3,1:3,iCamera)=U*V';
    dim(:,iCamera)=fscanf(fid,'%f',3);
    fclose(fid);
end
