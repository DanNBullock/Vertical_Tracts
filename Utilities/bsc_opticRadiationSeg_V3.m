function [RightMeyerFiber, RightMeyerBool, RightBaumFiber,RightBaumBool, LeftMeyerFiber, LeftMeyerBool, LeftBaumFiber,LeftBaumBool] =bsc_opticRadiationSeg_V3(wbfg, fsDir)
%[RightMdLF, RightMdLFindexes, LeftMdLF, LeftMdLFindexes] =bsc_segmentMdLF_neo(wbfg, fsDir)
%
% This function automatedly segments the middle longitudinal fasiculus
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
% -RightMdLF: fiber structure for right middle longitudinal fasiculus
% -LeftMdLF: fiber structure for left middle longitudinal fasiculus

% -RightMdLFIndexes: fiber indexes into the given wbfg for the right middle longitudinal fasiculus
% -LeftMdLFIndexes: fiber indexes into the given wbfg for the left middle longitudinal fasiculus

% (C) Daniel Bullock, 2017, Indiana University

%% parameter note & initialization

%maybe play with this if something isn't to your liking, it corresponds to
%the smoothing kernel used for the parietal and temporal ROIs
smoothParameter=5;
thalamusLut=[10 49];
choroidLut=[31 63];
blankBool(1:length(wbfg.fibers))=false;

    RightMeyerBool=blankBool';
    RightBaumBool=blankBool';
    LeftMeyerBool=blankBool';
    LeftBaumBool=blankBool';

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    

    
    %% thalamic ROI
    % we can be generous here, as we will be limiting our selection in
    % subsequent steps
    
 
    
    [thalamicROI] =wma_roiFromFSnums(fsDir,thalamusLut(leftright),1,7);
    thalamicROI.name='thalamicROI';
    
    anteriorThalamicPoint=max(thalamicROI.coords(:,2));
    lateralThalamicPoint=max(abs(thalamicROI.coords(:,1)))-3;
    postThalamicPoint=min(thalamicROI.coords(:,2));
    anteriorPoints=find(thalamicROI.coords(:,2)<anteriorThalamicPoint-0);
    %thalamicROI.coords=thalamicROI.coords(anteriorPoints ,: );
    
    
    
    %% posterior Occpital
    %generates the roi posterior occipital regions
    
    
    %119
    [postOccp] =bsc_roiFromFSnums(fsDir,[145 111 122 120]+sidenum,1,smoothParameter);
    postOccp.name='posteriorOccipitalROI';
    
    
    
    %% not ROIs
    % Make planar ROI 10 mm above the top voxel of the thalamic ROI
    topThalamicPoint=max(thalamicROI.coords(:,3));
    atlasNifti = wma_getAsegFile(fsDir , '2009');
        
    [topNotROI] =bsc_makePlanarROI(atlasNifti,topThalamicPoint+10,'z');
    topNotROI.name='topNotROI'; 
    
    % posterior bound of thalamus
    
    
    % prevent cross hemispheric fibers
    [midsaggitalNot] =bsc_makePlanarROI(atlasNifti,0,'x');
    midsaggitalNot.name='midsaggitalNot'; 
    
    %prevent anterior stretching fibers
    
    
    [antNotROI] =bsc_makePlanarROI(atlasNifti,anteriorThalamicPoint-5,'y');
    antNotROI.name='antNotROI'; 
    
    
    [precuneusROI] =bsc_roiFromFSnums(fsDir,[130]+sidenum,1,9);
 
    
    %% segmenting
    
    %set operands for ROIS
    operands={'endpoints','endpoints' 'not', 'not', 'not', 'not'};
    
    %switch for correct name
    if leftright==2
        sideflag='R';
    else
        sideflag='L';
    end
    currentFascicleName=strcat(sideflag,'_OR');
    
    %create object containing all rois
    currentROIs= [{postOccp} {thalamicROI} {topNotROI} {midsaggitalNot} {antNotROI}, {precuneusROI}];
    
    %actually segment
    [fascicle, FiberBoolVec] = wma_SegmentFascicleFromConnectome(wbfg, currentROIs, operands, currentFascicleName);
   
    %reporient streamlines
    fascicle = bsc_reorientFiber(fascicle);
    
    %determine if the 10th node is sufficiently far away from the first
    %node of the streamline in the y axis
    
    for iStreams=1:length(fascicle.fibers)
        
        %find node closest to posterior thalamic threshold
        postThalNode=find(min(abs(fascicle.fibers{iStreams}(2,:)-postThalamicPoint))==abs(fascicle.fibers{iStreams}(2,:)-postThalamicPoint));
        %check to see which side the fiber is on (fibers from other side
        %sometimes make it in due to inflated thalamic ROI
        if leftright==1
            expectedSign=-1;
        else
            expectedSign=1;
        end
        if expectedSign*fascicle.fibers{iStreams}(1,postThalNode)>0
            correctSideBool(iStreams)=true;
        else
            correctSideBool(iStreams)=false;
        end
        
        %using expected sign becaues of earlier use of abs
        if abs(fascicle.fibers{iStreams}(1,postThalNode)) >lateralThalamicPoint
        displacementVecBool(iStreams)=true;
        else
        displacementVecBool(iStreams)=false;
        end
    end
   
    
    MeyerFascicle=fascicle;
    BaumFascicle=fascicle;
    MeyerFascicle.fibers=fascicle.fibers(displacementVecBool&correctSideBool);
    BaumFascicle.fibers=fascicle.fibers(~displacementVecBool&correctSideBool);
    % BONUS SECRET:  THERE'S A TRACT CONNECTING THE THALAMUS TO ROI
    % 11/12130 OR MAYBE THE ANTERIOR CALCERINE SLUCUS
    % DOUBLE BONUS:  TO SEPARATE THE SUPERIOR AND INFERIOR BAUM LOOP, USE
    % THE HIPPOCAMPUS
    %% testing
    
%     figure
%     bsc_quickPlot(MeyerFascicle)
%     figure
%     bsc_quickPlot(BaumFascicle)

    
    %obtain fiber indexes corresponding to optik
    FiberIndexes=find(FiberBoolVec);
    MeyerIND=FiberIndexes(displacementVecBool);
    BaumIND=FiberIndexes(~displacementVecBool);
    clear displacementVecBool   
      
    %directs segmentation output to correct function output holder

    
    
    if leftright == 2
        RightMeyerFiber=MeyerFascicle;
        RightBaumFiber=BaumFascicle;
        
        RightMeyerBool(MeyerIND)=true;
        RightBaumBool(BaumIND)=true;
        
        RightBaumFiber.name='RBaum';
        RightMeyerFiber.name='RMeyer';

    else
        LeftMeyerFiber=MeyerFascicle;
        LeftBaumFiber=BaumFascicle;
        
        LeftMeyerBool(MeyerIND)=true;
        LeftBaumBool(BaumIND)=true;
   
        
        
            LeftBaumFiber.name='LBaum';
         LeftMeyerFiber.name='LMeyer';
    end
    
    clear displacementVecBool
    clear fascicle
    clear correctSideBool
end


  
    
    
    
    

end