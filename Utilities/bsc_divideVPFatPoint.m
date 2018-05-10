function [fgLat,fgMed, LatBool, MedBool] = bsc_divideVPFatPoint(fg, xCoordinate, yCoordinate)


%for playwithCoords=yCoordinate-9:3:yCoordinate+9
%yCoordinate=playwithCoords;
% fg=Left
% xCoordinate=xCoordinateL
% yCoordinate=yCoordinateL


for iFibers = 1:length(fg.fibers)
   yDists=abs(yCoordinate -fg.fibers {iFibers}(3,:));
   yDistMin=min(yDists);
   nodeOfInterest=find(yDistMin==yDists);
   
   TopPoints(:,iFibers)=fg.fibers {iFibers}(:,nodeOfInterest(1));
end



% firstMedialIndex=find(max(TopPoints(1,:))==TopPoints(1,:));
% firstLateralIndex=find(min(TopPoints(1,:))==TopPoints(1,:));
% 
% 
% fgMed.fibers{1}=fg.fibers{firstMedialIndex(1)};
% fg.fibers(firstMedialIndex(1))=[];
% fgLat.fibers{1}=fg.fibers{firstLateralIndex(1)};
% fg.fibers(firstLateralIndex(1))=[];



%% Now we move to sequentially checking whether adding the fiber to the left or the right maximizes the increase of distance between the average x value of the med and lat fiber groups



largerThanXIndexes=find(TopPoints(1,:)>xCoordinate);
lessThanXIndexes=find(TopPoints(1,:)<xCoordinate);


fgClear=fg;
fgClear.fibers=[];
fgMed=fgClear;
fgLat=fgClear;

largerThanXBool=(TopPoints(1,:)>xCoordinate);
lessThanXBool=(TopPoints(1,:)<xCoordinate);

if xCoordinate<0
    MedBool=largerThanXBool;
    LatBool=lessThanXBool;
    fgLat.name='LeftLatVPF';
    fgMed.name='LeftMedVPF';
    for imedfibers=1:length(lessThanXIndexes)
        fgLat.fibers{imedfibers,1}=fg.fibers{lessThanXIndexes(imedfibers)};
    end
    
    for ilatfibers=1:length(largerThanXIndexes)
        fgMed.fibers{ilatfibers,1}=  fg.fibers{largerThanXIndexes(ilatfibers)}  ;
    end
    %switches fgs because it gets them mixed up for left side
    

else
    LatBool=largerThanXBool;
    MedBool=lessThanXBool;
    fgLat.name='RightLatVPF';
    fgMed.name='RightMedVPF';
    for imedfibers=1:length(lessThanXIndexes)
        fgMed.fibers{imedfibers,1}=fg.fibers{lessThanXIndexes(imedfibers)};
    end
    
    for ilatfibers=1:length(largerThanXIndexes)
        fgLat.fibers{ilatfibers,1}=fg.fibers{largerThanXIndexes(ilatfibers)};
    end
end


% if xCoordinate<0
%     fgLat.name='LeftMedVPF';
%     fgMed.name='LeftLatVPF';
% else
%     fgLat.name='RightLatVPF';
%     fgMed.name='RightMedVPF';
% end
% 
% if xCoordinate>0
%     
%     fgLat.name='RightLatVPF';
%     fgMed.name='RightMedVPF';
% else
%     fgLat.name='LeftMedVPF';
%     fgMed.name='LeftLatVPF';
% end

% bsc_quickPlot(fgLat);
% view(180,0)
% bsc_quickPlot(fgMed);
% view(180,0)



end
%end

