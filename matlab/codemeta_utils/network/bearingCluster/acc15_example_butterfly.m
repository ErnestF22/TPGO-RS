function acc15_example_butterfly
fileDir='../figures/tikz';
fileNameBase='exampleButterfly';

for dataset=1:2
    switch dataset
        case 1
            [A,x]=bearingCluster_generateTest('butterfly');
        case 2
            [A,x]=bearingCluster_generateTest('butterfly-bend');
    end
    E.original=adj2edges(A,'oriented');
    E.a15=[E.original;1 5];
    E.a17=[E.original; 1 7];

    typesAddition=fields(E);

    for iType=1:length(typesAddition)
        type=typesAddition{iType};

        u.(type)=bearingCluster_getBearingsScalesFromE(x,E.(type));
        C.(type)=grOrientedCycleBasis(E.(type))';
        M.(type)=bearingCluster_measurementMatrixFromC(C.(type),u.(type),'methodNormalization','none');
        L.(type)=null(M.(type));
        [Lp.(type),l.(type)]=cnormalize(L.(type)');
        Lp.(type)=Lp.(type)';

        NAdd=size(E.(type),1)-size(E.original,1);
        ML.(type)=M.(type)*blkdiag(L.original,diag(ones(NAdd,1)));

        MCycle.(type)=cnormalize(ML.(type)(end-1:end,:));

        if ~strcmpi(type,'original')
            fileName=fullfile(fileDir,[fileNameBase num2str(dataset) '_' type '_data.tex']);
            fid=fopen(fileName,'wt');
            %fprintf(fid,'$\\setlength{\\arraycolsep}{2pt}\\bmat{');
            %fprintf(fid,' %+.3f & %+.3f & %+.3f\\\\', MCycle.(type));
            %fprintf(fid,'}$');
            %fclose(fid);
            disp(MCycle.(type))
            fprintf('Output written to %s\n',fileName);
            
            %bearingClustering_plot(x,E.(type),[ones(1,10) 2])
            fileName=fullfile(fileDir,[fileNameBase num2str(dataset) '_' type '.tex']);
            %bearingCluster_tikz(x,E.(type),'fileId',fileName,...
            %    'membership',[ones(1,10) 2],'flagNormalizeScale',false)

            ECycle=[1 2; 2 3; 3 1];
            xCycle.(type)=bearingCluster_embed(ECycle,MCycle.(type));
            %bearingClustering_plot(xCycle.(type),ECycle)
            fileName=fullfile(fileDir,[fileNameBase num2str(dataset) '_' type '_equiv.tex']);
            %bearingCluster_tikz(xCycle.(type),ECycle,'fileId',fileName, ...
            %    'scale',1)
        end
    end
end