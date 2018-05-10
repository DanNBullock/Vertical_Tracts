function [xTics, scaleType] = bsc_setFigureTics(dataStructure,idealSpace,scalePreference)

if notDefined('idealSpace')
    idealSpace=6;
end


dataMin=min(min(dataStructure));
dataMax=max(max(dataStructure));


switch scalePreference
    
    case 'linear'
        %higher than this and you are doing something wonky
        validIncriments=[5 10 25 100 250 500 1000 2000 5000 10000 20000];
        scaleType='linear';
        
        dataRange=dataMax-dataMin;
        targetRanges=5*validIncriments;
        incrimentIndex= targetRanges==min((targetRanges(targetRanges>dataRange)));
        
        expDiv=ceil(log10(dataRange)-1);
        MinVal=(floor((dataMin/(10^expDiv))))*(10^expDiv);
        MaxVal=(ceil((dataMax/(10^expDiv))))*(10^expDiv);
        xTics=MinVal:validIncriments(incrimentIndex):MaxVal;
        if length(xTics)<5
            xTics=MinVal:validIncriments(incrimentIndex):(MinVal+validIncriments(incrimentIndex)*5);
        end
        
    case 'log2'
        
        dataLog2Min=log2(dataMin);
        dataLog2Floor=floor(dataLog2Min);
        dataLog2Max=log2(dataMax);
        dataLog2Ceil=ceil(dataLog2Max);
        
        dataLog2Range=dataLog2Ceil-dataLog2Floor;
        
        for iUnits=1:5
            Log2spacings(iUnits)=ceil(dataLog2Range/iUnits);
            log2Start(iUnits)=(floor(dataLog2Floor/iUnits))*iUnits;
            log2End(iUnits)=(ceil(dataLog2Ceil/iUnits))*iUnits;
            log2WhiteSpace(iUnits)=abs(log2Start(iUnits)-dataLog2Min)+abs(log2End(iUnits)-dataLog2Max);
            log2WhiteSpaceRatio(iUnits)=log2WhiteSpace(iUnits)/(log2End(iUnits)-log2Start(iUnits));
            log2_idealspacing(iUnits)=(abs(1-Log2spacings(iUnits)/idealSpace)+1)^2;
            log2Fit(iUnits)=log2_idealspacing(iUnits)*log2WhiteSpaceRatio(iUnits);
        end
        
        log2FitMinIndex=find(min(log2Fit)==log2Fit);
        
        scaleType='log2';
        ticPowers=log2Start(log2FitMinIndex):log2FitMinIndex:log2End(log2FitMinIndex);
        for ipowers= 1:length(ticPowers)
            xTics(ipowers)=2^ticPowers(ipowers);
        end
        
    case 'log10'
        
        dataLog10Min=log10(dataMin);
        dataLog10Floor=floor(dataLog10Min);
        dataLog10Max=log10(dataMax);
        dataLog10Ceil=ceil(dataLog10Max);
        
        dataLog10Range=dataLog10Ceil-dataLog10Floor;
        
        for iUnits=1:5
            Log10spacings(iUnits)=ceil(dataLog10Range/iUnits);
            log10Start(iUnits)=(floor(dataLog10Floor/iUnits))*iUnits;
            log10End(iUnits)=(ceil(dataLog10Ceil/iUnits))*iUnits;
            log10WhiteSpace(iUnits)=abs(log10Start(iUnits)-dataLog10Min)+abs(log10End(iUnits)-dataLog10Max);
            log10WhiteSpaceRatio(iUnits)=log10WhiteSpace(iUnits)/(log10End(iUnits)-log10Start(iUnits));
            log10_idealspacing(iUnits)=(abs(1-Log10spacings(iUnits)/idealSpace)+1)^2;
            log10Fit(iUnits)=log10_idealspacing(iUnits)*log10WhiteSpaceRatio(iUnits);
        end
              
        log10FitMinIndex=find(min(log10Fit)==log10Fit);
        
        scaleType='log10';
        ticPowers=log10Start(log10FitMinIndex):log10FitMinIndex:log10End(log10FitMinIndex);
        for ipowers= 1:length(ticPowers)
            xTics(ipowers)=10^ticPowers(ipowers);
        end
end
end
