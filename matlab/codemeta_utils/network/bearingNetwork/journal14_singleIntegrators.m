function journal14_singleIntegrators
%generates data for plots in the journal paper
figDir='../../../papers/control/notes/figures';
figDim=[300,200];
flagSaveFigure=2;

allCostNamesBearings={'cosine'};%,'angleSq'};
costNameRanges='squared';

TShort=10;
TLong=1500;
TMotion=7500;

testSet='standard';

% 'b' = bearing only
% 'bp' = bearing+preconditioner
% 'bd' = bearing+distance
% 0 = no leaders
% 1 = one leader
% 2 = two leaders
switch testSet
    case 'standard'
        allConditions={{'b',0,TLong},{'bd',0,TLong}};
    case 'preconditioner'
        allConditions={{'b',0,TShort},{'bp',0,TShort}};
    case 'motion'
        allConditions={{'b',1,TMotion},{'bd',1,TMotion},{'b',2,TMotion}};
    case 'experimental'
        allConditions={{'bp',2,TMotion}};
end

t_cost.funsRanges=bearingCostFunctions(costNameRanges);

for d=2:3
    resetRands(1)
    switch d
        case 2
            t_node=bearingNetworkBuildTestNetwork();
            x0=6*randn(d,t_node.NNodes);
        case 3
            t_node=bearingNetworkBuildTestNetwork(11,3);
            x0=4*randn(d,t_node.NNodes);
    end
    
    offset=zeros(d,1);
    x0=x0-(mean(x0,2)+offset)*ones(1,t_node.NNodes);
    t_node.Ti=x0;
    t_node0=t_node;

    NNodes=t_node.NNodes;

    EBearings=t_node.E;
    ygBearings=t_node.Yijtruth;

    ERanges=t_node.Er;
    ygRanges=t_node.Yrijtruth;
    nygRanges=t_node.nYrijtruth;

    NCostBearings=length(allCostNamesBearings);
    NConditions=length(allConditions);

    save([mfilename '_' num2str(d) 'D_' testSet '_data'])

    for iCostBearings=1:NCostBearings
        costNameBearings=allCostNamesBearings{iCostBearings};
        disp(['# Cost ' costNameBearings])

        t_cost.funsBearings=bearingCostFunctions(costNameBearings);

        for iCondition=1:NConditions
            %make sure that each test uses the same initial conditions
            %(but with leaders the initial configuration is changed)
            t_node=t_node0;
            condition=allConditions{iCondition}{1};
            NLeaders=allConditions{iCondition}{2};
            TFinal=allConditions{iCondition}{3};
            baseFileName=[mfilename '_'];
            dxLeader=@(t) [-5/1000;0.01*sin(t/2/1000*pi)];
            if d==3
                baseFileName=[baseFileName '3D_'];
                dxLeader=@(t) [dxLeader(t);0];
            end
            baseFileName=[baseFileName costNameBearings '_' condition '_' ];
            if TFinal~=TLong
                baseFileName=[baseFileName 'T' num2str(TFinal) '_'];
            end
            t_cost.alpha=[1 1];
            switch condition
                case 'b'
                    disp('## Bearing formation')
                    t_cost.flagUseRanges=false;
                    optsControl={};
                case 'bp'
                    disp('## Bearing+preconditioner formation')
                    T=bearingNetworkPreconditioner(t_node,t_cost.funsBearings,'hop1n');
                    optsControl={'preconditioner',T,'preconditionerDelay',5};
                case 'bd'
                    disp('## Bearing+distance formation')
                    t_cost.flagUseRanges=true;
                    optsControl={};
                otherwise
                    error('Test condition for cost not recognized');
            end
            switch NLeaders
                case 0
                    %nothing to do
                case 1
                    idxLeader=4;
                case 2
                    idxLeader=[3 4];
                otherwise
                    error('Number of leaders invalid')
            end
            if NLeaders>0
                t_node.Ti=t_node.Titruth;
                optsControl=[optsControl {'leader',idxLeader,@(t) dxLeader(t)*ones(1,NLeaders)}];
                t_cost.alpha=[5 5];
                baseFileName=[baseFileName 'L' num2str(NLeaders) '_'];
            end            

            figure(1)
            %argSolver={'odeSolver',@odeEuler,'optsSolver',{'MaxStep',.1}};
            argSolver={'odeSolver',@ode15s};
            [t,x,t_node]=bearingNetworkEvolve(t_node,'tFinal',TFinal,'t_cost',t_cost,...
                'optsControl',optsControl,argSolver{:},...
                'showOdeProgress');
            output=bearingNetworkEvolveStats(t_node,t,x,'t_cost',t_cost,...
                'cost','angles','configurationDistance','residuals','rangedifference');

            save([baseFileName 'data'])

            %produce figures
            figBaseFileName=fullfile(figDir,baseFileName);
            figure(2)
            if NLeaders==0
                bearingNetworkPlot(t_node)
            else
                bearingNetworkPlot(t_node,'leader',idxLeader)
            end            
            hold on
            bearingNetworkPlotTrajectories(x)        
            hold off
            savefigure([figBaseFileName 'trajectories'],'epsc',figDim,flagSaveFigure)

            figure(3)
            semilogy(t,output.phi)
            savefigure([figBaseFileName 'cost'],'epsc',figDim,flagSaveFigure)

            figure(4)
            plot(t,output.m)
            savefigure([figBaseFileName 'centroid'],'epsc',figDim,flagSaveFigure)

            figure(5)
            plot(t,output.a)
            savefigure([figBaseFileName 'angles'],'epsc',figDim,flagSaveFigure)

            figure(6)
            semilogy(t,output.dist)
            savefigure([figBaseFileName 'configDistance'],'epsc',figDim,flagSaveFigure)

            figure(7)
            plot(t,output.c)
            savefigure([figBaseFileName 'resBearings'],'epsc',figDim,flagSaveFigure)

            if t_cost.flagUseRanges
                figure(8)
                plot(t,output.q)
                savefigure([figBaseFileName 'resRanges'],'epsc',figDim,flagSaveFigure)
            end
        end
    end
end