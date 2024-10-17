function packageSE3ExpLog

opts.packageName='SE3ExpLog';
opts.pairs={...
    'RiemannianGeometry/','./'...
    };
opts.baseDir='~/scratch/SE3ExpLog';
opts.autorightsPreambleFile='~/rvidal/cameranetw/generic/licensePreamble.txt';
opts.autorightsOpts={'dateVersion','authorName','Roberto Tron (tron@cis.jhu.edu)'};
opts.additionalFiles={'../utility/depPackage/LICENSE'};

depPackage('log_se3_test.m',opts)