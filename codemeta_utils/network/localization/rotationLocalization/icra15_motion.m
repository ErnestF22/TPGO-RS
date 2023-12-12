function icra15_motion
flagFigSave=2;
figDir='~/Documents/UPenn/papers/network/ICRA15-FrobeniousRotationLocalization/figures/';

load rotLocFrobConsensus_testMotion_data.mat
fs=setFigFontSize(8);
fn=setFigFont('Times');

legendText=[repmat('# Iterations = ',NNIter,1) num2str(NIter','%2d')];
eMean=rad2deg(squeeze(mean(eAligned)));
plot(t,eMean')
legend(legendText)
xlabel('s')
ylabel('Mean error [deg]')

disp('Steady-state errors')
disp(eMean(end,:))
savefigure(fullfile(figDir,'icra15_motion_errors'),'epsc',0.9*[270 180],flagFigSave);

setFigFontSize(fs);
setFigFont(fn);
