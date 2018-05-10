function [tractStats] = wma_singleTractAnalysis(indices,fe,WBFG,dt6)
%[summaryStructure] = wma_multiTractAnalysis(tractNameList,fiberIndexList,WBFG,FiberDir,saveHeader,dt6,fsDIR)
%
% OVERVIEW: this function runs a bevy of analysis functions including 
% bsc_multiVirtualLesion, bsc_multiTractVolumeFE, bsc_multiTractLengthFG,% (C) Daniel Bullock, 2017, Indiana University
% and bsc_multiTractProfilesFE.  It combines the output from several of
% these into a single summary structure and saves it down.
%
% INPUTS:
% -classification: Either the path to structure or the structure itself.
%  The strucure has a field "names" with (N) names of the tracts classified
%  while the field "indexes" has a j long vector (where  j = the nubmer of
%  streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%  a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%  indicatate that the streamline has been classified as a member of tract
%  (N).
%
% -WBFG: a wbFG object. 
%
% -fe: an fe object
%
% -fsDIR: path  to THIS SUBJECT'S freesurfer directory
%
% -dt6: either a path to a dt6 file or a dt6 object
%
% OUTPUTS:
% -summaryStructure: an amalgamated data structure containing the outputs
%  from bsc_multiVirtualLesion, bsc_multiTractLengthFG, and
%  bsc_multiTractProfilesFE.  See those respective functions for more
%  details.
%
% (C) Daniel Bullock, 2017, Indiana University
%% preliminaries

% conditional to prevent emptys from making it far
if ~isempty(indices) 

%% virtual lesion
%if it is a fe structure, get the wbFG out of it
if ~notDefined('fe')&&~isempty(fe)
        [vl] =wma_singleVirtualLesion(fe,indices);
        
        posIndicies=indices(ismember(indices,(find(fe.life.fit.weights>0))));

    else
        fprintf('\n fe structure not detected, will not run virtual lesion \n')

end
tractStats.vl=vl;

%% Quantative Stats
if ~isempty(posIndicies)
    [tractStats.morph.meanLength, tractStats.morph.lengthSTDEV, tractStats.morph.totalLength,...
        tractStats.morph.streamCount, tractStats.morph.volume,...
        tractStats.morph.streamLengths]=wma_singleTractStatQuantification(WBFG,posIndicies);
    tractStats.morph.lengthPerMM3=tractStats.morph.totalLength/tractStats.morph.volume;
    tractStats.morph.ValidProp=tractStats.morph.streamCount/length(posIndicies);
    
    %% Tract Profiles
    %doesn't work right now?
    if length(posIndicies)>5
        [tractStats.diffusion.fa, tractStats.diffusion.md, tractStats.diffusion.rd,...
            tractStats.diffusion.ad]=wma_singleTractProfile(posIndicies,WBFG,dt6);
    else
        tractStats.diffusion.fa=[];
        tractStats.diffusion.md=[];
        tractStats.diffusion.rd=[];
        tractStats.diffusion.ad=[];
    end
    
    tractStats.vl.evidencePerMM3=vl.evidence.s.mean/tractStats.morph.volume;
else
    tractStats=[] ;
end
else
    tractStats=[];
end
end