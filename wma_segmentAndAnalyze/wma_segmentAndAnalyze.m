function [results, classificationRAW]= wma_segmentAndAnalyze(fe,dt6,fsDIR)
%
% [results, classification]= wma_segmentAndAnalyze(fe,dt6,fsDIR)
%
% Overview:  This code takes in an fe structure (or FG) and then performs
% segmentation and analysis on both (the fe structure and the segmentation)
%
% Inputs:
%
% -fe: either a LiFE object or a path to one.
%
% -dt6: either a path to a dt6 file or a dt6 object
%
% -fsDIR: path  to THIS SUBJECT'S freesurfer directory
%
% Outputs:
%
% results: a summary structure containing information about both the
% connectome and the segmented tracts.
%%  preliminary loading
%
tic
if ischar(fe)
    load(fe);
    %haven't been able to test this with (path to wbFG passed in) yet.
    %if it is a fe structure, get the wbFG out of it
end

[classificationRAW]=wma_wrapper_PaperVersion(fe,dt6,fsDIR);


% differential pruning of AFQ and WMA tracts
classificationCut= removeOutliersClassification(classificationRAW,fe, 3, 3, 1:18);
classificationCut= removeOutliersClassification(classificationCut,fe, 4, 4, 21:34);

[tractStats] = wma_multiTractAnalysis(classificationCut,fe,dt6);

[~, results]= bsc_feAndAFQqualityCheck(fe, classificationCut);

if exist(strcat(fsDIR,'/wm_mask.nii.gz'),'file')
    fprintf('\n wm mask found')
    wmMasknii=niftiRead(strcat(fsDIR,'/wm_mask.nii.gz'));
    wmSize=sum(sum(sum(wmMasknii.data)));
    wmSize=wmSize*wmMasknii.pixdim(1)*wmMasknii.pixdim(2)*wmMasknii.pixdim(3);
    results.LiFEstats.WBFG.volume=wmSize;
    fprintf('\n wm volume is %i cubic mm', wmSize )
    
    for iFibers=1:length(tractStats)
        if ~isempty(tractStats{1,iFibers})
        tractStats{1,iFibers}.morph.volumeProportion=tractStats{1,iFibers}.morph.volume/results.LiFEstats.WBFG.volume;
        end
        %do something lke this for evidence as well.
    end
end

results.AFQstats.tractStats=tractStats;
classificationCut=wma_clearNonvalidClassifications(classificationCut,fe);
results.AFQstats.classification=classificationCut;
segTime=toc;
fprintf('\n Segmentation and analysis for %s has taken %4.2f hours.', fe.name, segTime/(60*60))
end

