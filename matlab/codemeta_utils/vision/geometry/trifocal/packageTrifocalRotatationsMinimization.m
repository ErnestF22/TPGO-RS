function packageTrifocalRotatationsMinimization

opts.packageName='trifocalMinimization';
opts.pairs={...
    'vision/geometry/','./';...
    };
opts.baseDir='~/scratch/trifocalMinimization';
opts.autorightsPreambleFile='~/Documents/JHU/vision/code/vision/cameraNetworks/utility/depPackage/licensePreamble.txt';
opts.autorightsOpts={'dateVersion','authorName','Roberto Tron (tron@seas.upenn.edu)'};
opts.additionalFiles={'~/Documents/JHU/vision/code/vision/cameraNetworks/utility/depPackage/LICENSE'};

depPackage({'~/Documents/UPenn/code/vision/geometry/POCmanoptTrifocalRotations.m'},opts)
