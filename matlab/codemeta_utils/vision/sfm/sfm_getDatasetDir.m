%Look up user name and return directory name for a named dataset
%function imgDirName=sfm_getDatasetDir(datasetName)
function imgDirName=sfm_getDatasetDir(datasetName)
user=getUser();
switch user
    case 'tron'
        datasetDirName='~/Documents/UPenn/datasets/';
        switch datasetName
            case {'fountain','castle','castleentry','herzjesu'}
                imgDirName=fullfile(datasetDirName, [datasetName '_dense']);
            case 'castlelarge'
                imgDirName=fullfile(datasetDirName, 'castle_dense_large');
            case 'herzjesularge'
                imgDirName=fullfile(datasetDirName, 'herzjesu_dense_large');
            case 'desk'
                imgDirName='~/Documents/UPenn/datasets/epipolarGeometry/desk/';
            otherwise
                error('datasetName not recognized')
        end
    case 'spiros'
        switch datasetName
          case 'fountain'
                imgDirName='~/Datasets/SFM/Fountain-P11';
            otherwise
                error('datasetName not recognized')
        end
    otherwise
        error('User %s not recognized. Please edit %s to add this user.',user,mfilename)
end
