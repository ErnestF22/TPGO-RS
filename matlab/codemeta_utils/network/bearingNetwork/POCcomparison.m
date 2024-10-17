function POCcomparison
resetRands()
TFinal=1000;
d=2;

methodBearing='gradient';
%methodBearing='projection';
%methodBearing='geodesic';
%methodBearing='franchi';
switch methodBearing
    case 'gradient'
        t_cost.funsBearings=bearingCostFunctions('cosine');
    case 'geodesic'
        t_cost.funsBearings.f=@(x) x^2/2;
    case {'projection','franchi'}
        t_cost.funsBearings=bearingCostFunctions('angle');
    otherwise
        error('methodBearing not recognized')
end

optsControl={'methodBearing',methodBearing};

t_node=bearingNetworkBuildTestNetwork();
if strcmp(methodBearing,'franchi')
    [t_node.Yijtruth,t_node.Gammaitruth,t_node.E]=bnFranchi_computeBearingAndGammaFormation(t_node.Titruth);
end
x0=6*randn(d,t_node.NNodes);
offset=zeros(d,1);
x0=x0-(mean(x0,2)+offset)*ones(1,t_node.NNodes);
t_node.Ti=x0;



argSolver={'odeSolver',@ode15s};
%argSolver={'odeSolver',@odeEuler};

figure(1)
switch methodBearing
    case {'gradient','geodesic','projection'}
        [t,x,t_node,ode]=bearingNetworkEvolve(t_node,'tFinal',TFinal,'t_cost',t_cost,...
            'optsControl',optsControl,argSolver{:},...
            'showOdeProgress');
    case 'franchi'
        [t,x,t_node,ode]=bnFranchi_Evolve(t_node,'tFinal',TFinal,argSolver{:},...
            'showOdeProgress');
end
        
output=bearingNetworkEvolveStats(t_node,t,x,'t_cost',t_cost,...
    'configurationDistance','angles');
outputControl=bearingNetworkEvolveStatsControl(ode,t,x,'totalNorm');

figure(2)
cla
bearingNetworkPlot(t_node)
hold on
bearingNetworkPlotTrajectories(x)        
hold off

save([mfilename '_data'])
