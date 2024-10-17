function acc15_example_cycles
fileDir='../figures/tikz';
fileNameBase='exampleCycle';

tikzOpts={'flagNormalizeScale',false};
for c=1:4
    switch c
        case 1
            A=adjgallery(3,'kneigh',1);
            E=adj2edges(A,'oriented');
            x=[ 0 0;
                3 0;
                2 1.5]';

            %bearingClustering_plot(x,E)
            fileName=fullfile(fileDir,[fileNameBase num2str(c) '.tex']);
            bearingCluster_tikz(x,E,'fileId',fileName,tikzOpts{:},'membership',[2 1 1])

        case 2
            A=adjgallery(4,'kneigh',1);
            E=adj2edges(A,'oriented');
            x=[ 0 0;
                3 0;
                2 1.5;
                0 1]';
            [L,u]=xEToLu(x,E);

            l2=bearingCluster_scalesSolutionFromL(L,-0.5,1,3);
            x2=bearingCluster_scales2nodes(E,u,l2);
            l3=bearingCluster_scalesSolutionFromL(L,0,1,3);
            x3=bearingCluster_scales2nodes(E,u,l3);
            %bearingClustering_plot(x,E)
            %hold on
            %bearingClustering_plot(x3,E)
            %hold off
            fileName=fullfile(fileDir,[fileNameBase num2str(c) 'a.tex']);
            bearingCluster_tikz(x,E,'fileId',fileName,tikzOpts{:},'membership',[2 1 1 1])
            fileName=fullfile(fileDir,[fileNameBase num2str(c) 'b.tex']);
            bearingCluster_tikz(x2,E,'fileId',fileName,tikzOpts{:},...
                'formatNameNode','vb','formatEdge','shadow','formatNode','vertex shadow')
            fileName=fullfile(fileDir,[fileNameBase num2str(c) 'c.tex']);
            bearingCluster_tikz(x3,E,'fileId',fileName,tikzOpts{:},...
                'formatNameNode','vb','formatEdge','shadow','formatNode','vertex shadow')

        case 3
            A=adjgallery(4,'kneigh',1);
            E=adj2edges(A,'oriented');
            x=[ 0 0 1;
                3 0 1;
                2 1.5 0.5;
                0 1 1]';

            %disp(size(xEToLu(x,E),2))
            %bearingClustering_plot(x,E)

            fileName=fullfile(fileDir,[fileNameBase num2str(c) '.tex']);
            bearingCluster_tikz(x,E,'fileId',fileName,tikzOpts{:},'membership',[2 1 1 1])
        case 4
            A=adjgallery(5,'kneigh',1);
            E=adj2edges(A,'oriented');
            x=[ 0 0 1;
                3 0 1;
                2 1.5 0.5;
                0 1 1;
                0 1 0]';
            [L,u]=xEToLu(x,E);
            disp(size(L,2))
            l2=bearingCluster_scalesSolutionFromL(L,1,1,3);
            x2=bearingCluster_scales2nodes(E,u,l2)+x(:,1)*ones(1,size(x,2));
            l3=bearingCluster_scalesSolutionFromL(L,0.5,1,3);
            x3=bearingCluster_scales2nodes(E,u,l3)+x(:,1)*ones(1,size(x,2));

            %bearingClustering_plot(x,E)
            %hold on
            %bearingClustering_plot(x2,E)
            %bearingClustering_plot(x3,E)
            %hold off
            fileName=fullfile(fileDir,[fileNameBase num2str(c) 'a.tex']);
            bearingCluster_tikz(x,E,'fileId',fileName,tikzOpts{:},'membership',[2 1 1 1 1])
            fileName=fullfile(fileDir,[fileNameBase num2str(c) 'b.tex']);
            bearingCluster_tikz(x2,E,'fileId',fileName,tikzOpts{:},...
                'formatNameNode','vb','formatEdge','shadow','formatNode','vertex shadow')
            fileName=fullfile(fileDir,[fileNameBase num2str(c) 'c.tex']);
            bearingCluster_tikz(x3,E,'fileId',fileName,tikzOpts{:},...
                'formatNameNode','vb','formatEdge','shadow','formatNode','vertex shadow')
            

    end
end


function [L,u]=xEToLu(x,E)
u=bearingCluster_getBearingsScalesFromE(x,E);
M=bearingCluster_measurementMatrixFromE(E,u);
L=null(M);
