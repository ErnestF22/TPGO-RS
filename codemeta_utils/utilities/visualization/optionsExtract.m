%Extracts only the named options from a cell array of optional parameters
%function optsExtracted=optionsExtract(opts,optsNames,optsNArgs,optsNamesExclude,optsNargsExclude)
%Inputs
%   opts                Cell array of options to parse (typically, opts)
%   optsNamesAndArgs    Cell array of options to extract. It must be a
%                       sequence of pairs 'OptionName',NumberOfArgs
%Optional inputs
%   optsNamesAndArgsExclude     Similar to optsNamesAndArgs, but these are
%                               specifically excluded
%NOTE: the parser will get confused if the argument to an unspecified
%option is the same as a specificed option name. To correct this situation, add
%the name of the unspecified option to optsNamesAndArgsExclude.

function optsExtracted=optionsExtract(opts,optsNamesAndArgs,optsNamesAndArgsExclude)
optsNames=optsNamesAndArgs(1:2:end);
optsNArgs=optsNamesAndArgs(2:2:end);
if exist('optsNamesAndArgsExclude','var')
    optsNamesExclude=optsNamesAndArgsExclude(1:2:end);
    optsNArgsExclude=optsNamesAndArgsExclude(2:2:end);
else
    optsNamesExclude={};
    optsNArgsExclude={};
end
optsExtracted={};
iOpts=1;
while(iOpts<=length(opts))
    flagJOpts=strcmpi(optsNames,opts{iOpts});
    if any(flagJOpts)
        nArgs=optsNArgs{flagJOpts};
        optsExtracted=[optsExtracted opts{iOpts:iOpts+nArgs}];              %#ok<AGROW>
        iOpts=iOpts+nArgs;
    end
    flagJOpts=strcmpi(optsNamesExclude,opts{iOpts});
    if any(flagJOpts)
        nArgs=optsNArgsExclude{flagJOpts};
        iOpts=iOpts+nArgs;
    end
    iOpts=iOpts+1;
end
