function bearingCostFunctions_test
argsList={{'power','cosineHalf',3},...
    {'cosine'},{'angle'},{'cosinesq'},...
    {'anglesq'},{'cosinel1l2'},{'anglel1l2'}};

for iArgs=5%1:length(argsList)
    funs=bearingCostFunctions(argsList{iArgs}{:});
    t=linspace(-0.9,1,201);
    figure
    subplot(3,1,1)
    check_der(funs.f,funs.df,t);
    title(argsList{iArgs}{1})
    subplot(3,1,2)
    check_der(funs.df,funs.ddf,t);
    subplot(3,1,3)
    plotfun(@(c) funs.f(c)+(1-c)*funs.df(c),t)
    xlabel('This needs to be negative for the convergence proof')
end