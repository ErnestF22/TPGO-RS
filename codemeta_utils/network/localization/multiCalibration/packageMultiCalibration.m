function packageMultiCalibration

opts.packageName='multiCalibration';
opts.pairs={...
    'network/localization/multiCalibration/','./';...
    'network/localization/anisotropic/','./';...
    };
opts.baseDir='~/scratch/multiCalibration';
opts.autorightsPreambleFile='~/Documents/JHU/vision/code/vision/cameraNetworks/utility/depPackage/licensePreamble.txt';
opts.autorightsOpts={'dateVersion','authorName','Roberto Tron (tron@seas.upenn.edu)'};
opts.additionalFiles={'~/Documents/JHU/vision/code/vision/cameraNetworks/utility/depPackage/LICENSE'};

depPackage({'~/Documents/UPenn/code/network/localization/multiCalibration/multiCalibrationNew_testRealData.m',...
    '~/Documents/UPenn/code/network/localization/multiCalibration/multiCalibration_test.m'},opts)
