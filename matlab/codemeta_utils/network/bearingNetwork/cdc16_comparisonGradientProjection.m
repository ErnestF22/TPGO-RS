function cdc16_comparisonGradientProjection
resetRands()

TLong=1000;
TMotion=7500;
t_cost.alpha=[5 5];

allTestSets={'standard','oneleader','twoleaders','twoleadersstar'}

for iTest=1:length(allTestSets)
    testSet=allTestSets{iTest};
    TFinal=TLong;
    flagMovingLeaders=false;
    switch testSet
        case 'standard'
            NLeaders=0;
            flagDoubleStarGraph=false;
            allMethodBearings={'gradient','gradientDistance','projection'};
        case 'oneleader'
            NLeaders=1;
            flagDoubleStarGraph=false;
            allMethodBearings={'gradient','gradientDistance','projection'};
        case 'twoleaders'
            NLeaders=2;
            flagDoubleStarGraph=false;
            allMethodBearings={'gradient','gradientDistance','projection'};
        case 'twoleadersstar'
            NLeaders=2;
            flagDoubleStarGraph=true;
            allMethodBearings={'franchi','gradientDistance','gradient','projection'};
    end

    optsControl={};
    t_cost.funsRanges=bearingCostFunctions('squared');

    for d=2:3
        resetRands(1)
        switch d
            case 2
                NNodes=7;
            case 3
                NNodes=11;
        end

        if ~flagDoubleStarGraph
            optsGraph={};
        else
            optsGraph={'Edges',bnFranchi_graphEdges(NNodes)};
        end
        t_node=bearingNetworkBuildTestNetwork(NNodes,d,optsGraph{:},'edgesRanges',[2 3; 3 2]);
        switch d
            case 2
                x0=6*randn(d,t_node.NNodes);
            case 3
                x0=4*randn(d,t_node.NNodes);
        end
        offset=zeros(d,1);
        x0=x0-(mean(x0,2)+offset)*ones(1,t_node.NNodes);


        NMethodBearings=length(allMethodBearings);

        t_cost.funsBearings=bearingCostFunctions('cosine');

        baseFileName=[mfilename '_'];
        switch d
            case 2
                baseFileName=[baseFileName '2D_'];
            case 3
                baseFileName=[baseFileName '3D_'];
        end

        switch NLeaders
            case 0
                %nothing to do
            case 1
                idxLeader=1;
            case 2
                idxLeader=[1 2];
            otherwise
                error('Number of leaders invalid')
        end
        if NLeaders>0
            x0(:,idxLeader)=t_node.Titruth(:,idxLeader);
            baseFileName=[baseFileName 'L' num2str(NLeaders)];
            if flagMovingLeaders
                dxLeader=@(t) [-5/1000;0.01*sin(t/2/1000*pi)];
                basfileName=[baseFileName 'M_'];
            else
                dxLeader=@(t) [0;0];
                basfileName=[baseFileName 'S_'];
            end
            if d==3
                dxLeader=@(t) [dxLeader(t);0];
            end            
            optsControl=[optsControl {'leader',idxLeader,@(t) dxLeader(t)*ones(1,NLeaders)}];
        end
        t_node.Ti=x0;
        t_node0=t_node;


        testSetFileName=[mfilename '_' num2str(d) 'D_' testSet];
        disp(['# Test set ' testSetFileName])

        save([testSetFileName '_data'])

        for iMethodBearings=1:NMethodBearings
            %make sure that each test uses the same initial conditions
            %(but with leaders the initial configuration is changed)
            t_node=t_node0;

            methodBearings=allMethodBearings{iMethodBearings};
            disp(['## Method: ' methodBearings])

            switch methodBearings
                case 'gradient'
                    t_cost.flagUseRanges=false;
                    optsControl=[optsControl {'methodBearing','gradient'}];
                case 'gradientDistance'
                    t_cost.flagUseRanges=true;
                    optsControl=[optsControl {'methodBearing','gradient'}];
                case 'projection'
                    t_cost.flagUseRanges=false;
                    optsControl=[optsControl {'methodBearing','projection'}];
                case 'franchi'
                    t_cost.flagUseRanges=false;
                    [t_node.Yijtruth,t_node.Gammaitruth,t_node.E]=bnFranchi_computeBearingAndGammaFormation(t_node.Titruth);
                otherwise
                    error('methodBearings not valid')
            end

            figure(1)
            %argSolver={'odeSolver',@odeEuler,'optsSolver',{'MaxStep',.1}};
            argSolver={'odeSolver',@ode15s};
            %argSolver={'odeSolver',@ode23tb};
            
            switch methodBearings
                case 'franchi'
                    [t,x,t_node,ode]=bnFranchi_Evolve(t_node,'tFinal',TFinal,argSolver{:},...
                        'showOdeProgress');
                otherwise
                    [t,x,t_node,ode]=bearingNetworkEvolve(t_node,'tFinal',TFinal,'t_cost',t_cost,...
                        'optsControl',optsControl,argSolver{:},...
                        'showOdeProgress');
            end
            results.output.(methodBearings)=bearingNetworkEvolveStats(t_node,t,x,'t_cost',t_cost,...
                'cost','angles','configurationDistance','residuals','rangedifference');
            results.outputControl.(methodBearings)=bearingNetworkEvolveStatsControl(ode,t,x,'totalNorm','cumulativeNormSq');
            results.t.(methodBearings)=t;
            results.x.(methodBearings)=x;
            results.t_node.(methodBearings).t_node=t_node;
        end

        save([testSetFileName '_data'])
    end
end