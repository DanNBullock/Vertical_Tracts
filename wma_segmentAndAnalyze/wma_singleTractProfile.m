function [fa, md, rd, ad]=wma_singleTractProfile(indices,wbFG,dt6)
% function [DiffusionProfile]=bsc_multiTractProfilesFE(tractNameList,fiberIndexList,wbFG,FiberDir,saveHeader,dt6, nosave)
%
% OVERVIEW:  This function runs dtiComputeDiffusionPropertiesAlongFG for
% the indexed streamlines from the input wbFG.  Doesn't load objects for
% expediancy's sake.
%
% INPUTS:
% -indices: a vector of indexes corresponding to the streamlines of
% interest.
%
% -wbFG:  a wbFG object.
%
% -FiberDir: directory path for the directory you would like your fiber
% indexes saved down to.
%
% -dt6: a dt6 object
%

%
% -OUTPUTS
% -DiffusionProfile: a structure wherein each field is named after a fiber tract.
% Each of these fields contains information like fractional anisotropy,
% mean diffusivity, radial diffusivity, and axial diffusivity.  If there 
% are no indexed fibers in % the corresponding fiberIndexList this will be
% blank (i.e. if segmentation failed).
%
% % (C) Daniel Bullock 2017 Bloomington
%
% maybe too short to justify its existance

%% preliminaries

% define tract
tract=wbFG;
tract.fibers=wbFG.fibers(indices);
    [fa, md, rd, ad, ~, ~, ~, ~, ~, ~]=dtiComputeDiffusionPropertiesAlongFG(tract, dt6, [] ,[] , 100);

end
