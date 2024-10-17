%An example that packages depPackage itself
function depPackage_test
%name of the package
opts.packageName='depPackage';
%define directory remappings
opts.pairs={...
    'depPackage/','./'...
    };
%staging directory relative to where the current script is
opts.baseDir=fullfile(mfilepath(),'..','..','releases','depPackage');
opts.autorightsOpts={'dateVersion','authorName','Roberto Tron (tron@bu.edu)','autoTag'};
%call depPackage with file(s) from which dependences should be computed
depPackage('depPackage_test.m',opts)
