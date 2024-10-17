function pdInit(name,N,flagOverwrite)
if ~exist('flagOverwrite','var')
    flagOverwrite=false;
end

indicatorsName=[name '_pdIndicators.mat'];

if ~exist(indicatorsName,'file') || flagOverwrite
    served=zeros(1,N);
    completed=zeros(1,N);
    save(indicatorsName, 'served', 'completed');
end