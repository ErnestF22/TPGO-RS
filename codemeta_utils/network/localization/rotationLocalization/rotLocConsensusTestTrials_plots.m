function rotLocConsensusTestTrials_plots
flagFigSave=2;
figDir='../../../../papers/network/IROS15-FrobeniousRotationLocalization/figures';

fs=setFigFontSize(8);
fn=setFigFont('Times');


dFrobenious=[];
dRiemannian=[];

%load('rotLocConsensusTestTrials_01-Mar-2015_14_54_15')
load('rotLocConsensusTestTrials_01-Mar-2015_16_29_10')

for iConfig=1:size(errors,1)
    for iTrial=1:size(errors,2)
        e=errors{iConfig,iTrial};
        if ~isempty(e)
            dFrobenious=[dFrobenious e.d.Frobenious];
            dRiemannian=[dRiemannian e.d.Riemannian];
        end
    end
end

semilogy([mean(dFrobenious,2) mean(dRiemannian,2)])
hold on
semilogy([median(dFrobenious,2) median(dRiemannian,2)],'--')
hold off
xlabel('Iterations')
ylabel('Aligned distance')
legend('Proposed','Intrinsic')

savefigure(fullfile(figDir,'convergence'),'epsc',0.9*[270 180],flagFigSave);
ax=axis();
ax(2)=5000;
axis(ax)

setFigFontSize(fs);
setFigFont(fn);
