function bsc_multiPlotData_v3(meanDataMatrix,stdDataMatrix,fiberNames,dataType,SaveDir,scalePref,groupmap)

%function bsc_multiPlotData_v3(meanDataMatrix,stdDataMatrix,fiberNames,dataType,SaveDir,scalePref,groupmap)
%
% CAN BE USED WITH APPROPRIATELY DESIGNED WRAPPER TO GENERATE MULTIPLE
% (INDIVIDUAL) PLOTS IN ONE GO
%
%   INPUTS
%
% > meanDataMatrix: A 3d matrix of mean data with subjects as rows, 
%                   replications by column, and tracts by series
%
% > stdDataMatrix: the same convention as meanDataMatrix except entries
%                  correspond to standard deviation for those measures
%
% > fiberNames: the tract names which will, in essence, serve as labels for
%               the series(s)
%
% > dataType: the type of data that is being plot, serves as label.
%
% > SaveDir:  directory you would like the figure saved down to.
%
% > scalePref:  prefered scaling for the x axis of the plot.  Accepts
%               either 'linear' 'log2' or 'log10'.  If nothing is entered,
%               will attempt to autodetect.
%
% > groupmap: numeric coding vector to indicate subjects' group membership.
%             i.e [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]
%
%   OUTPUTS
%
% [NONE], saves plot to disk.



%meanVec,stdVec,results.AFQstats.classification.names(tractsOfInterest),meanFields{iplots},SaveDir,[],groupmap



matrixDim=size(meanDataMatrix);
subjectNum=matrixDim(1);
fiberNum=matrixDim(2);

groups=unique(groupmap);

%if you need more you have a problem
%colorMapNames={'parula','autumn', 'winter', 'summer','spring','hot','cool'}
  
c{8} = colormap(parula(128));
c{1} = colormap(autumn(128));
c{7} = colormap(winter(128));
c{6} = colormap(summer(128));
c{5} = colormap(spring(128));
c{4} = colormap(hot(128));
c{2} = colormap(cool(128));
c{3} = colormap(hsv(128));

close all
colors=[];
for igroups=1:length(groups)
groupSubjects=sum(groupmap==igroups);
colorIndexes=1:floor((length(c{igroups})*(13/16))/groupSubjects):length(c{igroups})*(13/16);
colors=vertcat(colors,c{igroups}(colorIndexes,:));
end

span=.7;
delta = 0.8/(subjectNum-1);
incriments=0:span:span*subjectNum;
jitter = [incriments-(span*subjectNum/2)]*delta;

[xTics,logFlag] = bsc_getFigureTics_V5(meanDataMatrix);

if ~notDefined('scalePref')
    [xTics,logFlag]=bsc_setFigureTics(meanDataMatrix,[],scalePref);
end

boolVec(2:2:length(fiberNames))=true;

labelNames=fiberNames(boolVec);

close all
figure
hold on

%previously .5, 2 25
for isubjects= 1:subjectNum
    for Ifibers=1:fiberNum
        
        x =  Ifibers+ jitter(isubjects);
        
        switch logFlag
            case 'log10'
                
                semilogx(meanDataMatrix(isubjects,Ifibers),x,'o','markerfacecolor',colors(isubjects,:),  'markeredgecolor','k','linewidth',1,'markersize',35)
                
                semilogx([meanDataMatrix(isubjects,Ifibers) - stdDataMatrix(isubjects,Ifibers) ; ...
                    meanDataMatrix(isubjects,Ifibers) + stdDataMatrix(isubjects,Ifibers) ],  [x; x], '-','color',[.5 .5 .5],'linewidth',3)
            case 'linear'
                
                plot(meanDataMatrix(isubjects,Ifibers),x,'o','markerfacecolor',colors(isubjects,:),  'markeredgecolor','k','linewidth',1,'markersize',35)
                
                plot([meanDataMatrix(isubjects,Ifibers) - stdDataMatrix(isubjects,Ifibers) ; ...
                    meanDataMatrix(isubjects,Ifibers) + stdDataMatrix(isubjects,Ifibers) ],  [x; x], '-','color',[.5 .5 .5],'linewidth',3)
            case 'log2'
                
                plot(log2(meanDataMatrix(isubjects,Ifibers)),x,'o','markerfacecolor',colors(isubjects,:),  'markeredgecolor','k','linewidth',1,'markersize',35)
                
                plot([log2(meanDataMatrix(isubjects,Ifibers) - (stdDataMatrix(isubjects,Ifibers))) ; ...
                    log2(meanDataMatrix(isubjects,Ifibers) + (stdDataMatrix(isubjects,Ifibers))) ],  [x; x], '-','color',[.5 .5 .5],'linewidth',3)
        end
        
        
    end
    
    
end
fh=gcf;
fh.Name = strcat(dataType ,' - Major human whith matter tracts');
set(fh,'Position',[0,0,700,1050]);

switch logFlag
    case 'log10'
        set(gca, ...
            'xlim',[xTics(1) xTics(end)], 'xtick',(xTics), ...
            'ylim',[0 length(fiberNames)+1], 'ytick',[1.5:2:length(fiberNames)], ...
            'yticklabel',labelNames, ...
            'tickdir','out', ...
            'box','off', 'xscale','log',...
            'ticklen',[0.01 .01])
    case 'linear'
        
        set(gca, ...
            'ylim',[0 length(fiberNames)+1], 'ytick',[1.5:2:length(fiberNames)], ...
            'yticklabel',labelNames, ...
            'tickdir','out', ...
            'box','off', ...
            'ticklen',[0.01 .01])
    case 'log2'
        
        set(gca, ...
            'xlim',[log2(xTics(1)) log2(xTics(end))], 'xtick',log2((xTics)), ...
            'xticklabel',xTics, ...
            'ylim',[0 length(fiberNames)+1], 'ytick',[1.5:2:length(fiberNames)], ...
            'yticklabel',labelNames, ...
            'tickdir','out', ...
            'box','off', ...
            'ticklen',[0.01 .01])
end


figureNamed=(strcat(SaveDir,'_', dataType,'_Plot2.eps'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 30])
print ('-depsc', figureNamed ,'-r100')

end
