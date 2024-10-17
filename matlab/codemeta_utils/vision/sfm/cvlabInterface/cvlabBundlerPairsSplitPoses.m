function cvlabBundlerPairsSplitPoses(dirName)
listFiles=cvlabGetListFiles(dirName);
dirPairs=fullfile(dirName,'pairsKeysOnly');
s=load(fullfile(dirPairs,'pairsList'));
pairsList=s.pairsList;

for iPair=1:size(pairsList,1)
    nextPair=pairsList(iPair,:);
    dirNextPair=fullfile(dirPairs, [dirName 'pair_' num2str(nextPair(1)) '_' num2str(nextPair(2))]);
    
    dirNextPairCameras=fullfile(dirNextPair, 'cameras');
    disp(['Make dir ' dirNextPairCameras])
    
    fileName1=listFiles{nextPair(1)+1};
    fileName2=listFiles{nextPair(2)+1};
    
    srcFileName1=GetFullPath(fullfile(dirName,'cameras',fileName1));
    srcFileName2=GetFullPath(fullfile(dirName,'cameras',fileName2));
    
    cmd=['mkdir ' dirNextPairCameras];
    system(cmd);
    cmd=sprintf('ln -s %s %s/%s', srcFileName1, dirNextPairCameras, fileName1);
    system(cmd);
    cmd=sprintf('ln -s %s %s/%s', srcFileName2, dirNextPairCameras, fileName2);
    system(cmd);
end
