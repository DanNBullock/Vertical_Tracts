function  [L_pArc, L_TPC, L_pArc_Indexes,L_TPC_Indexes ,R_pArc, R_TPC, R_pArc_Indexes,R_TPC_Indexes] = bsc_automated_roi_segment_script_neo(wbfg, fsDir)
% [L_pArc, L_TPC, L_pArc_Indexes,L_TPC_Indexes ,R_pArc, R_TPC, R_pArc_Indexes,R_TPC_Indexes] = bsc_automated_roi_segment_script_neo(wbfg, fsDir)
%
% This function automatedly segments the posterior arcuate and temporo-parietal
% connection from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.
%
% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory
%
% Outputs:
% -L_pArc:fiber structure for left posterior arcuate
% -R_pArc:fiber structure for right posterior arcuate
% -L_TPC:fiber structure for left temporo-parietal connection
% -R_TPC:fiber structure for right temporo-parietal connection

% -L_pArc_Indexes:fiber indexes into the given wbfg for the left posterior arcuate
% -R_pArc_Indexes:fiber indexes into the given wbfg for the right posterior arcuate
% -L_TPC_Indexes:fiber indexes into the given wbfg for the left temporo-parietal connection
% -R_TPC_Indexes:fiber indexes into the given wbfg for the right temporo-parietal connection

% (C) Daniel Bullock, 2017, Indiana University


%% parameter note & initialization

    %this version relies on the aparc.a2009s+aseg file.  Here we make a
    %nii.gz version of one doesn't already exist.  keyboards out if there's
    %an error doing this
labelNifti= wma_getAsegFile(fsDir , '2009');
    
        %these 3 digit numbers correspond to the last 3 digits of the DK 2009
    %freesurfer look up table numbers.
    %(see:  https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT)
    %when added with sidenum you get the actual roi ID nums
    
    parietalROIs3=[157, 127, 168, 136, 126, 125];
    temporalROIs3=[121, 161, 137, 162, 138, 173];
    NOTROIs3=[122, 145, 166];

%maybe play with this if something isn't to your liking, it corresponds to
%the smoothing kernel used for the parietal and temporal ROIs

smoothParameter=5;

%% figure out a way to do this better or more elegantly
% these are basically hard set to just give us planes.  I could probably
% come up with a function to make planar ROIS if I wanted to...  for now
% this means we are no longer dependant on some saved down planar rois.
% Furthermore, the z=5, z=25 and y=-30 criteria could probably be
% customized to the particular subject, if that's what was desired.

% %formerly
% % pathToPlane1='/N/dc2/projects/lifebid/HCP/Dan/105115/dt6_b2000trilin/ROIs/Z_5_Plane.mat';
% % a plane ensuring that our vertical fibers are, in fact, running
% % vertically, maybe redundant given our current use of fs ROIS
% blankMatrix(1:255,1:255,1:255)=false;
% blankMatrix(53:203,13:213,133)=true;
% p1Indexes=find(blankMatrix);
% [p1X, p1Y, p1Z]=ind2sub(size(blankMatrix),p1Indexes);
% p1X=p1X-128; p1Y=p1Y-128; p1Z=p1Z-128;
% Plane1=dtiNewRoi('Z_5_Plane','r',[p1X, p1Y, p1Z]);

[Plane1]=bsc_makePlanarROI(labelNifti,5, 'z');

% %formerly
% % pathToPlane2='/N/dc2/projects/lifebid/HCP/Dan/105115/dt6_b2000trilin/ROIs/Z_25_Plane.mat';
% % a plane ensuring that our vertical fibers are, in fact, running
% % vertically, maybe redundant given our current use of fs ROIS
% blankMatrix(1:255,1:255,1:255)=false;
% blankMatrix(53:203,13:213,158)=true;
% p1Indexes=find(blankMatrix);
% [p1X, p1Y, p1Z]=ind2sub(size(blankMatrix),p1Indexes);
% p1X=p1X-128; p1Y=p1Y-128; p1Z=p1Z-128;
% Plane2=dtiNewRoi('Z_25_Plane','r',[p1X, p1Y, p1Z]);

[Plane2]=bsc_makePlanarROI(labelNifti,25, 'z');

%formerly
% pathToPlane3='/N/dc2/projects/lifebid/HCP/Dan/105115/dt6_b2000trilin/ROIs/Y_-30_Plane.mat';
% a plane ensuring that our candidate pArc or TPC fibers do not extend
% too far forward.  Some temporal ROIs extend quite anterior, so this
% criterion is necessary
% blankMatrix(1:255,1:255,1:255)=false;
% blankMatrix(53:203,98,68:213)=true;
% p1Indexes=find(blankMatrix);
% [p1X, p1Y, p1Z]=ind2sub(size(blankMatrix),p1Indexes);
% p1X=p1X-128; p1Y=p1Y-128; p1Z=p1Z-128;
% Plane3=dtiNewRoi('Y_-30_Plane','r',[p1X, p1Y, p1Z]);

[Plane3]=bsc_makePlanarROI(labelNifti,-25, 'y');

%formerly
% pathToPlane5='/N/dc2/projects/lifebid/HCP/Dan/105115/dt6_b2000trilin/ROIs/not_saggital_0X.mat';
% a plane ensuring that we don't get any fibers that are crossing the
% hemispehres, maybe not necessary but doesn't hurt anything.
% blankMatrix(1:255,1:255,1:255)=false;
% blankMatrix(128,13:213,68:213)=true;
% p1Indexes=find(blankMatrix);
% [p1X, p1Y, p1Z]=ind2sub(size(blankMatrix),p1Indexes);
% p1X=p1X-128; p1Y=p1Y-128; p1Z=p1Z-128;
% Plane4=dtiNewRoi('X_0_Plane','r',[p1X, p1Y, p1Z]);

[Plane4]=bsc_makePlanarROI(labelNifti,0, 'x');


%% documentation from previous versions
% pathToPlane4='/N/dc2/projects/lifebid/HCP/Dan/105115/dt6_b2000trilin/ROIs/NOT_box_ventricles_Z15.mat';
% this capacity was fufilled by the later generation of the NOT roi using
% fs indexes 4 and 43

% %pathToPlane5='/N/dc2/projects/lifebid/HCP/Dan/105115/dt6_b2000trilin/ROIs/not_sphere_0_-75_10_30.mat';
% this was typically used to prevent fibers from traversing around the
% medial side of the ventricles.  Turns out this isn't necessary because
% neither the TPC nor the pArc actually travel medial to the ventricles.
% Furthermore, I am unfamiliar with any reports of fiber tracts traveling
% vertically and medial to the ventricles.

%% actual segmentation

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
     
    %% parietal ROI
    %generates the roi for the parietal regions corresponding to the pArc
    %and TPC
    
    [mergedParietalROI] =bsc_roiFromFSnums(fsDir,parietalROIs3+sidenum,1,smoothParameter);
    mergedParietalROI.name='parietalROI';  
        
    %% temporal ROI
    %generates the roi for the temporal regions corresponding to the pArc
    %and TPC
    [mergedTemporalROI] =bsc_roiFromFSnums(fsDir,temporalROIs3+sidenum,1,smoothParameter);
    mergedParietalROI.name='temporalROI';
    
    %% NOT ROI
    %if things look weird, maybe smooth the ventricle ROI?
    [mergedNOTROI] =bsc_roiFromFSnums(fsDir,NOTROIs3+sidenum);
    mergedmergedNOTROIROI.name='mergedNOTROI';  

    %NOTE THAT WE DO NOT SMOOTH THE NOT ROI
    
    fprintf('\n ROIs genereated for side %i', leftright)
    
    %% optional debugging code for visualizing rois
    %     %testing/debugging purposes, plot to check where your rois are
    %     parietalIndexes=find(parietalROImask);
    %     [testx,testy,testz]=ind2sub(size(parietalROImask),parietalIndexes);
    %     figure
    %     scatter3(testx,testy,testz, 'MarkerEdgeColor','r')
    %     view([-90, 90])
    %% segment
    
    %set operands for ROIS
    operands={'endpoints','endpoints', 'and', 'and', 'not', 'not', 'not' };
    
    %switch for correct name
    if leftright == 2
        currentFascicleName=('right_parietal_TO_temporalFibers_pArcTPC_Candidate');
    else
        currentFascicleName=('left_parietal_TO_temporalFibers_pArcTPC_Candidate');
    end
    
    %create object containing all rois
    currentROIs= [{mergedParietalROI} {mergedTemporalROI} {Plane1} {Plane2} {Plane3} {Plane4} {mergedNOTROI}];
    
    %actually segment
    [fascicle, FiberBoolVec] = wma_SegmentFascicleFromConnectome(wbfg, currentROIs, operands, currentFascicleName);
    
    %obtain fiber indexes corresponding to the pArc + TPC amalgum
    FiberIndexes=find(FiberBoolVec);
    fprintf('\n %s segmented',currentFascicleName )
    
    %optional cleaning step.  Left out for now 
    if ~isempty(fascicle.fibers)
        %optional cleaning, left out here, typically we save down messy and
        % clean after reloading, right before an operation is performed
        % (i.e. plotting or virtual lesioning)
        %         [~, keep] = mbaComputeFibersOutliers(fascicle,4,4);
        %         fprintf('\n[%s] Found a tract with %i fibers... \n',mfilename,sum(keep))
        %         fprintf(' %s \n', currentFascicleName);
        %         fascicle = fgExtract(fascicle,find(keep),'keep');
    else
        fprintf('No fibers found for %s, segmentation error probable \n', currentFascicleName);
    end
    
    fascicle.name = currentFascicleName;
    
    
    %% Now split them
    % because the previous step returns an amagum of the TPC and pArc we
    % now have to divide them using the bottom point of the IPS as
    % indicated by freesurfer
    % 157 =S_intrapariet_and_P_trans_label
    
    %generate an ROI from the aparc aseg file.
    [bottom_of_IPS_ROI] =bsc_roiFromFSnums(fsDir,sidenum+157,1,smoothParameter);
    
    %this is a method of obtaining the appropriate z and x coordinates of
    %the bottom of the IPS.  We dont really care about the y coordinate for
    %the purposes of splitting the pArc and TPC.
    Zthresh=min(bottom_of_IPS_ROI.coords(:,3))+15;
    Zmatchindexes=find((bottom_of_IPS_ROI.coords(:,3)-(Zthresh-15))<.75);
    Xrange=unique(bottom_of_IPS_ROI.coords(Zmatchindexes,1));
    midX=mean(Xrange);
    uniqueY=unique(bottom_of_IPS_ROI.coords(:,2));
    midY=median(uniqueY);
    midlineY=find(bottom_of_IPS_ROI.coords(:,2)-midY<.75);
    midYZvals=unique(bottom_of_IPS_ROI.coords(midlineY,3));
    midZ=median(midYZvals);
    
    %function which divides the amalgum given the coordinates computed
    %above.  The VPF component of the function name is a holdover from
    %older versions of the code.
    [fgLat,fgMed, LatBool, MedBool] = bsc_divideVPFatPoint(fascicle, midX, midZ);
    
    %THE ABOVE CODE ACTUALLY DISTINGUISHES BETWEEN LEFT AND RIGHT, SEE IF
    %YOU CAN INCORPORATE THAT HERE.
    if leftright == 2
        R_pArc_Indexes=FiberIndexes(LatBool);
        R_TPC_Indexes=FiberIndexes(MedBool);
        R_pArc=fgLat;
        R_TPC=fgMed;
    else
        L_pArc_Indexes=FiberIndexes(LatBool);
        L_TPC_Indexes=FiberIndexes(MedBool);
        L_pArc=fgLat;
        L_TPC=fgMed;
    end
    
    
end
end
