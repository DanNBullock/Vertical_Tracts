function[figHandle, results]= bsc_feAndAFQqualityCheck(fe, classification, saveDir)
%function[figHandle, results]= bsc_feAndAFQqualityCheck(fe, afq_classificiation, saveDir)
%
% This function computes a number of statistics for an input fe structure,
% and, if included, an afq segmentation.  if you include a saveDir variable
% it saves down a figure illustrating a number of these statistics
% (relative to length) as well as the results structure.  "Validated" or
% "post-life" refers to streamlines which have been validated by life and
% thus have nonzero weights in fe.life.fit.weights.  Correspondingly,
% "WBFG" or "pre-life" referrs to ALL streamlines from the input whole
% brain tractography.
%
% INPUTS:
% -fe:  either a path to an fe structure, or an fe structure itself.
%
% -classification: Either the path to structure or the structure itself.
%  The strucure has a field "names" with (N) names of the tracts classified
%  while the field "indexes" has a j long vector (where  j = the nubmer of
%  streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%  a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%  indicatate that the streamline has been classified as a member of tract
%  (N).
%
% -saveDir: path that you would like the results struture and plot saved
%       to.  If not defined, does not save the output.
%
% OUTPUTS:
%
% -figHandle: handle for the mutiplot figure.  Can be used to adjust figure
%       properties or to save figure.
%
% -results: a structure with multiple fields summarizing various properties
%       of your input FE structure and (if included) AFQ segmentation.  See
%       code for specifics of all fields and their interpretation
%
% % (C) Daniel Bullock, 2017, Indiana University
%% preliminaries

% loads AFQ classification file if it is a path
if ~notDefined('classification')
    if ischar(classification)
        load(classification)
        %catches classification if it was saved.  Can probably be made more
        %robust to other saving schemas
    end
end

if ischar(fe)
    load(fe);
end

% gets streamlines from the FE structure
wbFG=feGet(fe,'fibers acpc');
% NOTE: the acpc transform must be applied to the WBFG because the fibers
% in the fe.fg structure are in image space (which would result in oddities
% when computing lengths)

for inames=1:length(classification.names)
leftNames(inames)=~isempty(strfind(classification.names{inames},'Left'));
rightNames(inames)=~isempty(strfind(classification.names{inames},'Right'));
interHemiNames(inames)=~leftNames(inames) & ~rightNames(inames);
end



% gets positively weighted streamlines and their indexes
posIndexes=find(fe.life.fit.weights>0);

%a test for the
if isempty(find(isnan(fe.life.fit.weights),1))>0
    nanNum=length(find(isnan(fe.life.fit.weights)));
    fprintf('\n %i NaN values detected in fe.life.fit.weights', nanNum)
end

pos_fg=wbFG;
pos_fg.fibers=wbFG.fibers(posIndexes);



%% resultStrucStuff
% computes length  of positively weighted streamlines

for istreamlines=1:length(pos_fg.fibers)
    pos_fg_streamLengths(istreamlines)=sum(sqrt(sum(diff(pos_fg.fibers{istreamlines},1,2).^2)));
   
end


% Computes length  of all streamlines
for istreamlines=1:length(wbFG.fibers)
    wbFG_streamLengths(istreamlines)=sum(sqrt(sum(diff(wbFG.fibers{istreamlines},1,2).^2)));
    
end


results.LiFEstats.fe_name=fe.name;

% Computes the average and standard deviation of the voxelwise error. See
% feGet: 'voxrmses0norm' for more information.
% Concern: does this include non tractography occupied voxels?
allvoxelsErr= feGet(fe, 'voxrmses0norm');
results.LiFEstats.RMSE.voxel_average_=mean(allvoxelsErr(~isnan(allvoxelsErr)));
results.LiFEstats.RMSE.voxel_stdev=std(allvoxelsErr(~isnan(allvoxelsErr)));
results.LiFEstats.RMSE.voxel_count=length(allvoxelsErr);
% Gets the total RMSE for the tractography model.  Not equivalent to
% sum(allvoxelsErr(~isnan(allvoxelsErr))), hence the above concern.
results.LiFEstats.RMSE.WB=feGet(fe,'rmsetotal');
% Whole model error sum
results.LiFEstats.RMSE.WB_norm_total=sum(allvoxelsErr(~isnan(allvoxelsErr)));

% Computes total length of pre and post LiFE streamline connectome.
% Compares the values to get proportional reduction in total connectome
% streamline length.
results.LiFEstats.WBFG.length_total=sum(wbFG_streamLengths);
results.LiFEstats.validated.streams_length_total=sum(pos_fg_streamLengths);
results.LiFEstats.validated.WBFG_length_reduc_proportion=results.LiFEstats.validated.streams_length_total/results.LiFEstats.WBFG.length_total;

%Volume Measures
%results.LiFEstats.WBFG.WMvolume=length(wbFGVolVec);
%results.LiFEstats.validated.WMvolume=length(posVolVec);

% Computes total number of connectome streamlines and the number that are
% validated by LiFE.  Computes the survivorship proportion.  Standard for
% HCP with single paramater tractography is ~ 15 to 20 percent.  Ensembel
% approaches push this to ~30.
results.LiFEstats.validated.stream_count=length(posIndexes);
results.LiFEstats.WBFG.stream_count=length(wbFG.fibers);
results.LiFEstats.validated.WBFG_proportion=results.LiFEstats.validated.stream_count/results.LiFEstats.WBFG.stream_count;

% avgerage and standard deviation for streamline lengths
results.LiFEstats.validated.avg_stream_length=mean(pos_fg_streamLengths);
results.LiFEstats.WBFG.avg_stream_length=mean(wbFG_streamLengths);
results.LiFEstats.validated.stream_length_stdev=std(pos_fg_streamLengths);
results.LiFEstats.WBFG.stream_length_stdev=std(wbFG_streamLengths);

% maximum streamline length.  Minimum wouldn't be interesting as this would
% simply be whatever your tractograhpy generation algorithm's minimum
% streamline lenght parameter was (i.e., 10 mm)
results.LiFEstats.validated.stream_max_length=max(pos_fg_streamLengths);
results.LiFEstats.WBFG.stream_max_length=max(wbFG_streamLengths);



%% AFQ additions
% if a classification (path or object) was passed into the function, adds
% AFQ relevant stats to the results structure, otherwise moves on.
if ~notDefined('classification')
    % Vector containging indexes of all streamlines classified by AFQ
    AFQIndexVec=find(classification.index);
    
    if length(classification.index)==results.LiFEstats.validated.stream_count
        
        life2afqFlag=true;
        
        % If the number of items in the AFQ classification structure is equal
        % to the number of positively weighted fibers in the WBFG (extracted
        % from the FE struc) this suggests that the user performed AFQ on the
        % positively weighted subset of streamlines from the WBFG structure,
        % rather than running it on the WBFG that was initially used to create
        % the FE structure.
        fprintf ('\n Number of items in AFQ classification structure suggests that AFQ was performed on a validated ONLY connectome WBFG.')
        fprintf ('\n Proceeding with this assumption.')
        % in this case we cannot know how many AFQ classified streamlines there
        % were before LiFE was run.
        results.AFQstats.WBFG.classified_stream_count=[];
        results.AFQstats.WBFG.classified_stream_proportion=[];
        results.AFQstats.WBFG.classified_stream_avg_length=[];
        results.AFQstats.WBFG.classified_stream_length_std=[];
        
        % Counts number of nonzero (i.e., classified) entries in the
        % classification.index field.
        results.AFQstats.validated.classified_stream_count=length(find(classification.index));
        % Divides the AFQ classification count by the total WBFG streamline
        % count to get the proportion of streamlines classified by AFQ
        results.AFQstats.validated.classified_stream_proportion=results.AFQstats.validated.classified_stream_count/results.LiFEstats.WBFG.stream_count;
        % NOTE: these streamlines are both AFQ classified AND LiFE validated
        
        % in theory this is unnecessary, as this vector should just be
        % every streamline. i.e. each member of AFQIndexVec is a member of
        % PosIndexes
        survivorVec=AFQIndexVec(ismember(AFQIndexVec,posIndexes));
        
        results.AFQstats.validated.classified_stream_avg_length=mean(wbFG_streamLengths(AFQIndexVec));
        results.AFQstats.validated.classified_stream_length_stdev=std(wbFG_streamLengths(AFQIndexVec));
        
    elseif length(classification.index)==results.LiFEstats.WBFG.stream_count
        
         life2afqFlag=false;
        
        % If the number of items in the AFQ classification structure is NOT
        % equal to the number of positively weighted fibers in the WBFG we
        % assume that AFQ was run independent of LiFE
        
        % Counts number of WBFG streamlines classified by AFQ and computes
        % proportion reltaive to total WBFG streamline count.
        results.AFQstats.WBFG.classified_stream_count=length(AFQIndexVec);
        results.AFQstats.WBFG.classified_stream_proportion=results.AFQstats.WBFG.classified_stream_count/results.LiFEstats.WBFG.stream_count;
        
        % Compute mean and standard deviation of AFQ classified streamlines
        % in the WBFG
        results.AFQstats.WBFG.classified_stream_avg_length=mean(wbFG_streamLengths(AFQIndexVec));
        results.AFQstats.WBFG.classified_stream_length_stdev=std(wbFG_streamLengths(AFQIndexVec));
        
        % Vector containing streamline indexes that were classified by AFQ and
        % validated by LiFE
        survivorVec=AFQIndexVec(ismember(AFQIndexVec,posIndexes));
        % Counts number of streamlines that were classified by AFQ and also
        % validated by LiFE.
        results.AFQstats.validated.classified_stream_count=length(survivorVec);
        results.AFQstats.validated.classified_stream_proportion=results.AFQstats.validated.classified_stream_count/results.LiFEstats.WBFG.stream_count;
        
        % Compute mean and standard deviation of validated and AFQ classified
        % streamlines in the WBFG
        results.AFQstats.validated.classified_stream_avg_length=mean(wbFG_streamLengths(survivorVec));
        results.AFQstats.validated.classified_stream_length_stdev=std(wbFG_streamLengths(survivorVec));
        
    else
        
        fprintf('\n Number of streamlines classified by AFQ is neither equal to total number of streamlines in WBFG nor total number of validated streamlines. \n')
        fprintf('\n Length of classification.index field should be exactly equal to number of streamlines from input tractography structure\n');
        % usually running afq with interhemispheric split ends up
        % generating about 800 extra fibers, which completely screws up the
        % indexing scheme.
        fprintf('\n Likely cause is default interhemisphericsplit behavior of AFQ\n');
        keyboard
        
    end
    
end
%% figure setup and computation
figure

% Sets figure into appropriate aspect ratio.  First two entries in [2 2 25 5]
% are arbitrary
set(gcf,'Units','inches',...
    'Position',[2 2 25 10]);

% Compute the bin counts of the WBFG and validatedstreamline lengths,
% binned at 1mm.
[WBFGhist, ~]=histcounts(wbFG_streamLengths,(10:220));
[validHist, ~]=histcounts(pos_fg_streamLengths,(10:220));

%% AFQ CASE
if ~notDefined('classification')
    % Do the same for the AFQ classified tracts
    [WBFGclassHist, ~]=histcounts(wbFG_streamLengths(AFQIndexVec),(10:220));
    [validClassHist, ~]=histcounts(wbFG_streamLengths(survivorVec),(10:220));
    
    %% standard Life Plots
    % Plot comparing pre and post life streamline counts by length
    subplot(3,3,1)
    hold on
    plot ((11:220),WBFGhist,'b', 'LineWidth',1.25)
    plot ((11:220),validHist,'r', 'LineWidth',1.25)
    title('WBFG & Validated Stream Count Comparison')
    legend('WBFG','Validated')
    xlabel('Streamline Length (mm)')
    ylabel('Streamline Count')
    
    % Plot comparing normalized pre and post life streamline counts by length
    subplot(3,3,2)
    hold on
    plot ((11:220),WBFGhist/results.LiFEstats.WBFG.stream_count,'b', 'LineWidth',1.25)
    plot ((11:220),validHist/results.LiFEstats.validated.stream_count,'r', 'LineWidth',1.25)
    title('Normalized WBFG & Validated Stream Count Comparison')
    legend('WBFG','Validated')
    xlabel('Streamline Length (mm)')
    ylabel('Whole Brain proportion')
    
    % Plot illustrating the bias in the validated streamlines as assocaited
    % with length
    subplot(3,3,3)
    hold on
    plot ((11:220),((WBFGhist/results.LiFEstats.WBFG.stream_count)-(validHist/results.LiFEstats.validated.stream_count))*10000,'g', 'LineWidth',1.25)
    plot ((11:220),(zeros(1,210))*10000,'r', 'LineWidth',1.25)
    title('Validation Bias Relative to Streamline Length')
    ylim([-40,25])
    legend('WBFG ratio - Validated ratio','No Bias')
    xlabel('Streamline Length (mm)')
    ylabel('Validation bias (%)')
    
    %% Plots Specific to AFQ Output
    % Plot comparing pre and post life streamline counts by length
    subplot(3,3,5)
    hold on
    plot ((11:220),WBFGclassHist,'b', 'LineWidth',1.25)
    plot ((11:220),validClassHist,'r', 'LineWidth',1.25)
    title('AFQ Classified WBFG & Validated Stream Count Comparison')
    legend('WBFG, AFQ classified','Validated & AFQ classified')
    xlabel('Streamline Length (mm)')
    ylabel('Streamline Count')
    
    % Plot comparing normalized pre and post life streamline counts by length
    subplot(3,3,6)
    hold on
    
    %plot ((11:220),WBFGclassHist/results.AFQstats.WBFG.classified_stream_count,'b', 'LineWidth',1.25)
    %plot ((11:220),validClassHist/results.AFQstats.validated.classified_stream_count,'r', 'LineWidth',1.25)
    plot ((11:220),(WBFGclassHist/results.LiFEstats.WBFG.stream_count)*100,'b', 'LineWidth',1.25)
    plot ((11:220),(validClassHist/results.LiFEstats.validated.stream_count)*100,'r', 'LineWidth',1.25)
    
    title('Normalized, AFQ Classified WBFG & Validated Stream Count Comparison')
    legend('WBFG, AFQ classified','Validated & AFQ classified')
    xlabel('Streamline Length (mm)')
    ylabel('Proportion of Whole Brain Streamlines Classified (%)')
    
%     % Plot illustrating the bias in the validated streamlines as assocaited
%     % with length
%     if ~life2afqFlag
%         % Computation gets unweildy, so here we use an intermediary variable
%         computeVec=((WBFGclassHist/results.AFQstats.WBFG.classified_stream_count)-(validClassHist/results.AFQstats.validated.classified_stream_count))-...
%             ((WBFGhist/results.LiFEstats.WBFG.stream_count)- (validHist/results.LiFEstats.validated.stream_count));
%         
%         subplot(2,3,6)
%         hold on
%         plot ((11:220),((WBFGclassHist/results.AFQstats.WBFG.classified_stream_count)-(validClassHist/results.AFQstats.validated.classified_stream_count)),'g', 'LineWidth',1.25)
%         plot ((11:220),computeVec,'c', 'LineWidth',1.25)
%         title('Validation bias relative to streamline length and specificity to AFQ classification')
%         ylim([-.0175,.01])
%         legend('WBFG ratio - Validated ratio, AFQ classified - General','WBFG ratio - Validated ratio, AFQ classified - Specific')
%         xlabel('Streamline Length (mm)')
%         ylabel('VAlidation bias')
%     end
cumValid=zeros(1,length(validHist));
cumWBFG=cumValid;
for ilengths=5:length(validHist)
cumValid(ilengths)=(validHist(ilengths)/results.LiFEstats.validated.stream_count)+sum(cumValid(ilengths-1));
cumWBFG(ilengths)= (WBFGhist(ilengths)/results.LiFEstats.WBFG.stream_count)+sum(cumWBFG(ilengths-1));
end

    subplot(3,3,4)
    hold on
    plot ((11:220),cumWBFG,'b', 'LineWidth',1.25)
    plot ((11:220),cumValid,'r', 'LineWidth',1.25)
    title('Cumulative portion of fibers in connectome, by length')
    legend('WBFG','Validated')
    xlabel('Streamline Length (mm)')
    ylabel('Portion of tracts less than or equal to length')

    
    %computation for afq bar plot
        tractNames=classification.names;
    for itracts=1:length(tractNames)
        tractProportion(itracts)=length(find(classification.index==itracts))/length(classification.index);
    end
    
    % THIS IS A KUDGE DUE TO THE ISSUES CAUSED BY INTERHEMISPHERIC FIBERS
    % AND ODD NUMBERS OF TRACTS
    plotInput=zeros(2,(length(tractProportion)-sum(interHemiNames))/2);
    
    plotInput(1,1)=tractProportion(1);
    plotInput(2,1)=sum(tractProportion( 2:sum(interHemiNames)));
     
    plotInput(1,:)=tractProportion((sum(interHemiNames)+1):2:end);
    plotInput(2,:)=tractProportion((sum(interHemiNames)+2):2:end);
    labelNames{1}=tractNames{1};
    for ilabels=2:(length(tractNames)-(sum(interHemiNames)))/2
        curName=tractNames{ilabels*2};
        spaceindexes=strfind(curName,' ');
        if ~isempty(spaceindexes)
        labelNames{ilabels}=curName(spaceindexes(1)+1:end);
        else
            labelNames{ilabels}=curName;
        end
    end
    
    bottomPlot=subplot(3,3,[7,8,9])
    hold on
    bar((plotInput')*100)
    title('Proportion of connectome streamlines in tract')
    legend('Left','Right')
    xlabel('Tract')
    ylabel('% classificaiton input streamlines in tract (%)')
    ylim([-0 2])
    set(gca,'xtick',[1:1:length(labelNames)])
    set(gca,'XTickLabel',labelNames, 'FontSize',8,'FontName','Times')
    bottomPlot.XTickLabelRotation=-45;
    
    %textual otputs (save for later)
    
    %text(.5,1.75)
    
else
    %% standard Life Plots
    % Plot comparing pre and post life streamline counts by length
    subplot(1,3,1)
    hold on
    plot ((11:220),WBFGhist,'b', 'LineWidth',1.25)
    plot ((11:220),validHist,'r', 'LineWidth',1.25)
    title('WBFG & Validated Stream Count Comparison')
    legend('WBFG','Validated')
    xlabel('Streamline Length (mm)')
    ylabel('Streamline Count')
    
    % Plot comparing normalized pre and post life streamline counts by length
    subplot(1,3,2)
    hold on
    plot ((11:220),WBFGhist/results.LiFEstats.WBFG.stream_count,'b', 'LineWidth',1.25)
    plot ((11:220),validHist/results.LiFEstats.validated.stream_count,'r', 'LineWidth',1.25)
    title('Normalized WBFG & Validated Stream Count Comparison')
    legend('WBFG','Validated')
    xlabel('Streamline Length (mm)')
    ylabel('Whole Brain proportion')
    
    % Plot illustrating the bias in the validated streamlines as assocaited
    % with length
    subplot(1,3,3)
    plot ((11:220),((WBFGhist/results.LiFEstats.WBFG.stream_count)-(validHist/results.LiFEstats.validated.stream_count)),'g', 'LineWidth',1.25)
    title('Validation Bias Relative to Streamline Length')
    ylim([-.0175,.01])
    legend('WBFG ratio - Validated ratio')
    xlabel('Streamline Length (mm)')
    ylabel('Validation bias')
    
end
%gets figure handle
figHandle=gcf;

if ~notDefined('saveDir')
    saveas(gcf,strcat(saveDir,'/fe_AFQ_StatPlot.eps'));
    save(strcat(saveDir,'/fe_AFQ_resultStruc.mat'),'results')
end

end


