function [fg, keep]=  bsc_tractByEndpointROIs(fg, rois)


minDist = 0.87;




  
    %cut streamlines to just the first and last node
    for istreamlines=1:length(fg.fibers)
       
        endpoint1(:,istreamlines)=fg.fibers{istreamlines}(:,1) ;
        endpoint2(:,istreamlines)=fg.fibers{istreamlines}(:,end);

    end
    
    
    

        [~, distr1e1]=nearpoints(endpoint1, rois{1}.coords');
        [~, distr1e2]=nearpoints(endpoint2, rois{1}.coords');
        
        [~, distr2e1]=nearpoints(endpoint1, rois{2}.coords');
        [~, distr2e2]=nearpoints(endpoint2, rois{2}.coords');
        
        
        keep=or(and(distr1e1<minDist,distr2e2<minDist),and(distr2e1<minDist,distr1e2<minDist));
        
        
       
            
       fg.fibers=fg.fibers(keep);
end