function [vl] =wma_singleVirtualLesion(fe,indices)
% function [FullVL] =bsc_singleVirtualLesion(fe,indicies)
%
% OVERVIEW:  singleton variant of  bsc_multiVirtualLesion.  Made for quick
% and straightforward compuation of virtual lesioning.  Probably redundant
% with other code
%
% INPUTS:
%
% -fe:  the FE structure that the virtual lesioning will be conducted upon.
%  it is from the FE's tractome that the streamlines will be removed and
%  the effect computed.
%
% -indicies: A vector conatining integer indicies corresponding to what is
%  (presumed to be) a tract of interest.
%
% OUTPUTS:
% -vl: structure with information like evidence, headerData (from
%  saveHeader), tract name (from the tract list), and other data
%  corresponding to the virtual lesion.  If there are no indexed fibers in
%  the corresponding fiberIndexList this will be blank.
%
% % (C) Daniel Bullock 2017 Bloomington

%% Virtual lesion
if length(indices)>4
[ vl.rmse_wVL, vl.rmse_woVL, vl.nFib_tract, vl.nFib_PN, vl.nVoxels ] = feComputeVirtualLesion_norm(fe, indices);
% feVirtualLesion for plots
fevl = feComputeEvidence(vl.rmse_woVL, vl.rmse_wVL);
vl.evidence=fevl;
else
vl=[];
% why even bother with this?
end