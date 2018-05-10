function [meanLength, lengthSTDEV, totalLength, count, volume, streamLengths] =wma_singleTractStatQuantification(wbFG,indices)
% function [meanLength, lengthSTDEV, totalLength, count, volume, streamLengths] =wma_singleTractStatQuantification(wbFG,indices)
%
% OVERVIEW: a singleton variant of wma_multiTractStatQuantification.  For
% quick and straightforward computation of tract stats.
%
% INPUTS:
%
% -wbFG:  a whole brain fiber group structure. 
%
% -indicies: a vector of integers corresponding to the indexes of what is
%  (presumably) a tract that is being assessed quantatatively.
%
% OUTPUTS:
% -meanLength: average tract length
%
% -lengthSTDEV: standard deviation of length
%
% -totalLength: sum of all streamline lengths for tract
%
% -streamLengths=vector containing computed lengths of each streamline
%
% -count= total count of streamlines
%
% -volume: volume occupied by tract
%
% % (C) Daniel Bullock 2017 Bloomington, Indiana
%% calculation of statistics of interest
volVec=[];
tractFibers=wbFG.fibers(indices);
for istreamlines=1:length(tractFibers)
    streamLengths(istreamlines)=sum(sqrt(sum(diff(tractFibers{istreamlines},1,2).^2)));
    volVec=horzcat(volVec,tractFibers{istreamlines});
end
meanLength=mean(streamLengths);
lengthSTDEV=std(streamLengths);
count=length(streamLengths);
%why doesn't it like totalLength?
totalLength=sum(streamLengths);

% volume, as indicated by number of 1mm voxels occupied by at least 1 tract
% node.  In theory if you wanted to resample this computation (i.e. how
% many .5 mm voxels are occpupied by at least 1 tract node) you could
% just multiply the volVec by a factor.  I.e. .5mm -> 2 or .25 mm -> 4
volume=length(unique(floor(volVec'),'rows'));

end
