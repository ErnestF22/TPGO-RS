function sfm_utilityCopyDatasets
datasets={'fountain','castle','castleentry','castlelarge','herzjesu','herzjesularge'};
for d=datasets;
    fileName1=['sfm_test_data_' d{1} '.mat'];
    fileName2=['sfmdata_' d{1} '.mat'];
    disp(['Copying ' fileName1 ' to ' fileName2])
    copyfile(fileName1,fileName2)
end
