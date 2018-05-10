function [tractStats] = wma_multiTractAnalysis(classification,feORwbFG,dt6)
%[summaryStructure] = wma_multiTractAnalysis(tractNameList,fiberIndexList,feORwbFG,FiberDir,saveHeader,dt6,fsDIR)
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
% -feORwbFG: either a path to a saved wbFG or a wbFG object.  If you pass a fe
% structure or an fe path it will (probably) extract the fg from it
% appropriately and subsequently used the positively weighted tracts when relevant.
%
% -dt6: either a path to a dt6 file or a dt6 object
%
% -fsDIR: path  to THIS SUBJECT'S freesurfer directory
%
% OUTPUTS:
% tractStats: an amalgamated data structure containing the outputs
% from wma_singleVirtualLesion, wma_singleTractStatQuantification, and
% wma_singleTractProfile.  See those respective functions for more
% details.
%
% (C) Daniel Bullock, 2017, Indiana University
%% preliminaries

% parses fe/wbfg input
[WBFG, fe] = bsc_LoadAndParseFiberStructure(feORwbFG);


% loads the dt6 file if a path was passed
if ischar(dt6)
    dt6 = dtiLoadDt6(dt6);
end

%% iterative tract analysis


for Itracts=1:length(classification.names)
    indices=find(classification.index==Itracts);
    [tractStats{Itracts}] = wma_singleTractAnalysis(indices,fe,WBFG,dt6);
    
end
end