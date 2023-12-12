%Packages the algorithms for QuickShift and QuickMatch
function quickShiftMatching_package
%name of the package
opts.packageName='quickShiftMatching';
%define directory remappings
opts.pairs={...
    'learning/quickshift','./';...
    'learning/','./';...
    'math/utilities/','./utilities';...
    'math/EuclideanGeometry/','./utilities';...
    'math/','./'...
    };
%staging directory relative to where the current script is
opts.baseDir=fullfile(mfilepath(),'..','..','releases','quickShiftMatching');
opts.autorightsOpts={'dateVersion','authorName','Roberto Tron (tron@bu.edu)','autoTag'};
%additional files
opts.additionalFiles={fullfile(mfilepath(),'README.md')};

%call depPackage with file(s) from which dependences should be computed
depPackage({'quickshift_test.m','quickshift_matching_test.m'},opts)
