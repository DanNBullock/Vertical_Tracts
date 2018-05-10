function [L_VOF, R_VOF, L_VOF_Indexes, R_VOF_Indexes] = ...
       wma_segment_vof(fg, fsROIdir,classification, dt, thresh,v_crit, arcThresh, parcThresh)
% Segment the VOF from a wholebrain connectome
%
% [L_VOF, R_VOF, L_VOF_Indexes, R_VOF_Indexes] = wma_segment_vof(wholebrainfgPath,,fsROIdir,thresh,v_crit, dt, arcThresh, parcThresh)
%
% This function will take in a wholebrain connectome, a segmented arcuate
% fasciculus and a freesurfer segmentation and return the vertical
% occipital fasciculus (VOF).
%
% Inputs:
%
% wholebrainfgPath - A path (or fg structure) for a wholebrain fiber group.
% fsROIdir         - Path to a directory containing the aparcAseg file for
%                    THIS subject.  The standard FS subject directory
% thresh           - A fiber must travel vertical for a large proportion of
%                    its length. The default values are likely fine
% vcrit            - To find fibers that we can considder vertical, we must
%                    define a threshold of how far a fiber must travel in
%                    the vertical direction. v_crit defines how much
%                    farther a fiber must travel in the vertical direction
%                    compare to other directions (e.g., longitudinal) to be
%                    a candidate for the VOF. The default is 1.3
% dt               - dt6.mat structure or path to the dt6.mat file.
%
% Outputs
% L_VOF, R_VOF     - Left and right hemisphere VOF fiber groups

% Copyright Jason D. Yeatman, September 2014. Code released with:
% Yeatman J.D., Weiner K.S., Pestilli F., Rokem A., Mezer A., Wandell B.A.
% (2014). The vertical occipital fasciculus: A forgotten highway. PNAS.
%
%  Edited by Franco Pestilli and Daniel Bullock, 2017
%
%% Argument checking and parameter setting
% Path to ROIs

% Remove any fiber that doesn't go vertical (positive z) for thresh% of its
% coordinates
if notDefined('thresh')
    thresh = [.95 .6];
end
% Fibers must travel this much farther vertically than other directions
if notDefined('v_crit')
    v_crit = 1.3;
end

% Default is to define VOF as fibers that have fewer than 20 nodes of
% overlap with the arcuate
if notDefined('arcThresh')
    arcThresh = 20;
end
% Default is to consider fibers that are anterior to the posterior arcuate
% as not part of the VOF
if notDefined('parcThresh')
    parcThresh = 1;
end

%% Find vertical fibers

% this will almost certianly cause a problem if it is used elsewhere
% with a different naming schema.

for iTracts = 1:length(classification.names)
    if strcmp(classification.names{iTracts}(1),'L')
        leftVec(iTracts)=true;
        rightVec(iTracts)=false;
    else
        leftVec(iTracts)=false;
        rightVec(iTracts)=true;
    end
    if ~isempty(strfind(classification.names{iTracts},'Arcuate'))
        ArcuateVec(iTracts)=true;
    else
        ArcuateVec(iTracts)=false;
    end
    if ~isempty(strfind(classification.names{iTracts},'pArc'))
        pArcVec(iTracts)=true;
    else
        pArcVec(iTracts)=false;
    end
end

L_arcuate  = dtiNewFiberGroup('L_arcuate');
L_arcIndexes=find(classification.index==(find(leftVec&ArcuateVec)));
L_arcuate.fibers=fg.fibers(L_arcIndexes);

R_arcuate  = dtiNewFiberGroup('R_arcuate');
R_arcIndexes=find(classification.index==(find(rightVec&ArcuateVec)));
R_arcuate.fibers=fg.fibers(R_arcIndexes);

L_pArc  = dtiNewFiberGroup('L_pArc');
L_pArcIndexes=find(classification.index==(find(leftVec&pArcVec)));
L_pArc.fibers=fg.fibers(L_pArcIndexes);

R_pArc  = dtiNewFiberGroup('L_pArc');
R_pArcIndexes=find(classification.index==(find(leftVec&pArcVec)));
R_pArc.fibers=fg.fibers(R_pArcIndexes);

% From the wholebrain fiber group find all the vertical fibers that
% terminate in ventral occipitotemporal cortex (VOT).

[L_fg_vert, R_fg_vert, L_vertical_fascicles_identities, R_vertical_fascicles_identities] = wma_find_vertical_fibers(fg,fsROIdir);


%% Separate VOF from arcuate
L_VOF      = dtiNewFiberGroup('L_VOF');
L_pArc_vot = dtiNewFiberGroup('L_posteriorArcuate_vot');
R_VOF      = dtiNewFiberGroup('R_VOF');
R_pArc_vot = dtiNewFiberGroup('R_posteriorArcuate_vot');

if ~isempty(L_fg_vert.fibers)
    
    % Make an arcuate fiber density image
    arcFdImg = dtiComputeFiberDensityNoGUI(L_arcuate,dt.xformToAcpc,size(dt.b0));
    % Theshold image at voxels with >n fibers
    arcFdImg = single(arcFdImg>1);
    % Compute the number of nodes that overlap with the arcuate
    fgVals = dtiGetValFromFibers(arcFdImg,L_fg_vert,inv(dt.xformToAcpc));
    fgMvals = cellfun(@(x) sum(x),fgVals);
    % Remove fibers that overlap with the arcuate for more than arcThresh nodes
    L_VOF.fibers = L_fg_vert.fibers(fgMvals<arcThresh);
    L_VOF_Indexes = L_vertical_fascicles_identities(fgMvals<arcThresh);

    
    % From the VOF fiber group, remove any fibers that are further anterior
    % than the core of the posterior arcuate.
    if parcThresh == 1
        ymaxL      = mean(cellfun(@(x) mean(x(2,:)),L_pArc.fibers));
        L_VOF_keep = cellfun(@(x) all(x(2,:) < ymaxL),L_VOF.fibers);
        
        % Add fibers that are being removed to pArc_vot

        L_VOF.fibers  = L_VOF.fibers(L_VOF_keep);
        L_VOF_Indexes = L_VOF_Indexes(L_VOF_keep);
    end
    
else
    L_VOF = [];

end

%% Repeat for the right hemisphere
if ~isempty(R_fg_vert.fibers)
    
    % Make an arcuate fiber density image
    arcFdImg = dtiComputeFiberDensityNoGUI(R_arcuate,dt.xformToAcpc,size(dt.b0));
    % Theshold image at voxels with >2 fibers
    arcFdImg = single(arcFdImg>2);
    
    % Compute the number of nodes that overlap with the arcuate
    fgVals  = dtiGetValFromFibers(arcFdImg,R_fg_vert,inv(dt.xformToAcpc));
    fgMvals = cellfun(@(x) sum(x),fgVals);
    
    % Remove fibers that overlap with the arcuate for more than 20 nodes
    R_VOF.fibers       = R_fg_vert.fibers(fgMvals<arcThresh);
    R_VOF_Indexes      = R_vertical_fascicles_identities(fgMvals<arcThresh);
    
    
    % From the VOF fiber group, remove any fibers that are further anterior
    % than the core of the posterior arcuate.
    if parcThresh == 1
        ymaxR      = mean(cellfun(@(x) mean(x(2,:)),R_pArc.fibers));
        R_VOF_keep = cellfun(@(x) all(x(2,:)<ymaxR),R_VOF.fibers);
        
        % Add fibers that are being removed to pArc_vot

        R_VOF.fibers       = R_VOF.fibers(R_VOF_keep);
        R_VOF_Indexes      = R_VOF_Indexes(R_VOF_keep);
    end
    
else
    R_VOF=[];

end

return