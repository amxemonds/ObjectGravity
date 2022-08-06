
close all;
clear all;

load('horizon0_catmullData.mat')
optDegDiff = 12;

%% 0 HORIZON

numCells = numCells8;
recordAccumulationCells = recordAccumulationCells8;

for cell = 1:1:numCells

    if cell < 9
        range2 = [-50,-37.5,-25,-12.5,0,12.5,25,37.5,50];   % (50 to -50)
        whichComparisonIndices = [1,26,51,76,101]; % e.g. 0, 12.5, 25, 37.5, 50
        rangeObjOri = 1:9;
        catromRef = 62.5;
        context = 'original';

    elseif cell <= 43
        range2 = [-100,-75,-50,-25,0,25,50,75,100];         % (100 to -100)
        whichComparisonIndices = [1,51,101,151,201,251,301]; % e.g. -50, -25, 0, 25, 50, 75, 100
        rangeObjOri = 2:10;
        catromRef = 125;
        context = 'revised';
    else
        range2 = [-50,-37.5,-25,-12.5,0,12.5,25,37.5,50];   % (50 to -50)
        whichComparisonIndices = [1,26,51,76,101]; % e.g. 0, 12.5, 25, 37.5, 50
        rangeObjOri = 1:9;
        catromRef = 62.5;
        context = 'original';
    end
   
    % these are labeled -50 to 50, but are actually 50 to -50

    % have one distribution for MT25 : 8p is recordAccumulationCells8 2
    meanMT25raw = recordAccumulationCells(:,1,cell,2); % may either have 9 or 11 values

    % have one distribution for MT-25 : 8n is recordAccumulationCells8 1
    meanMTn25raw = recordAccumulationCells(:,1,cell,1); % may either have 9 or 11 values

    if (sum(meanMT25raw)==0) || (sum(meanMTn25raw)==0)
        continue;
    end

    countedCell = 1;

    % whichProcess,'smoothBoxcar': meanMT25fit and meanMTn25fit have 9 values
    % find the values at the smooth fit for MT25
    meanMT25rawPadded = padarray(meanMT25raw(rangeObjOri),[1,0],'replicate');
    meanMT25fit = conv(meanMT25rawPadded',[1/3,1/3,1/3],'valid')';

    % find the values at the smooth fit for MT-25
    meanMTn25rawPadded = padarray(meanMTn25raw(rangeObjOri),[1,0],'replicate');
    meanMTn25fit = conv(meanMTn25rawPadded',[1/3,1/3,1/3],'valid')';

    % catmull-rom splines for shifted analysis
    rangePadded = [-catromRef,range2,catromRef]';

    meanMT25fitPadded = padarray(meanMT25fit,[1,0],'replicate');
    C25 = createCatmullRomSpline([rangePadded,meanMT25fitPadded],context);
    C25 = C25(2:end-1,:);

    for nn = 2:length(range2)-1
        idx = find(round(C25(:,1),1)==range2(nn));

        if idx
            C25(idx(1),:) = []; 
        end
    end

    C25 = unique(C25,'rows');
    
    meanMTn25fitPadded = padarray(meanMTn25fit,[1,0],'replicate');
    Cn25 = createCatmullRomSpline([rangePadded,meanMTn25fitPadded],context);
    Cn25 = Cn25(2:end-1,:);

    for nn = 2:length(range2)-1
        idx = find(round(Cn25(:,1),1)==range2(nn));

        if idx
            Cn25(idx(1),:) = []; 
        end
    end

    Cn25 = unique(Cn25,'rows');
    
    step = optDegDiff+1;

    if contains(context,'original')
        % retinal hypothesis
        % original
        % n: -50 to 0
        % p: 0 to +50

        % 1:101, 101:201
        % 2:102, 100:200
        % 3:103, 99:199

        % ret
        n25 = step:(length(C25)-1)/2+step;
        p25 = (length(C25)-1)/2+2-step:length(C25)+1-step;
        
        % grav
        p25Default = (length(C25)-1)/2+2-1:length(C25)+1-1;
        n25Default = (length(C25)-1)/2+2-1:length(C25)+1-1;

    elseif contains(context,'revised')
        % retinal hypothesis
        % revised
        % n: -50 to 100 
        % p: -100 to 50

        % 1:301, 101:401
        % 2:302, 100:400
        % 3:303, 99:399

        % ret
        n25 = step:(length(C25)-1)/2+100+step;
        p25 = (length(C25)-1)/2-100+2-step:length(C25)+1-step;
        
        % grav
        p25Default = (length(C25)-1)/2-100+2-1:length(C25)+1-1;
        n25Default = (length(C25)-1)/2-100+2-1:length(C25)+1-1;
    end
    
    gcf = figure;
    allPlotted25 = [meanMT25raw(rangeObjOri);meanMT25fit];
    allPlottedn25 = [meanMTn25raw(rangeObjOri);meanMTn25fit];

    % raw data, MT25
    retX = C25(p25,1);
    retX = retX(whichComparisonIndices);

    retY = C25(p25,2);
    retY = retY(whichComparisonIndices);

    gravX = C25(p25Default,1);
    gravX = gravX(whichComparisonIndices);

    gravY = C25(p25Default,2);
    gravY = gravY(whichComparisonIndices);

    % boxcar smoothed data, MT25
    subplot(2,1,2)
    hold on;

    yyaxis left
    barHandles = length(range2);

    for rangeIdx = 1:1:length(range2)
        barHandles(rangeIdx) = bar(range2(rangeIdx),meanMT25fit(rangeIdx),'BarWidth', 5,'LineWidth',0.1,'EdgeColor','w');
        set(barHandles(rangeIdx), 'FaceColor', [0.9,0.9,0.9],'EdgeColor','none');
    end

    axis([range2(1)-10,range2(end)+10,min(allPlotted25)-2,max(allPlotted25)+2])
    xlabel('Object Lean in Degrees','fontname','lato','FontSize', 15)
    ylabel('Firing Rate','fontname','lato','FontSize', 15)
    title('+25º Monkey Tilt, Boxcar Smoothed Data, 0º Horizon','fontname','lato','FontSize',15)

    yyaxis right
    plot(C25(:,1),C25(:,2),'LineWidth',2,'color','k');
    plot(gravX,gravY,'v','LineWidth',2,'color','m','MarkerSize',20,'MarkerFaceColor','m');
    plot(retX,retY,'^','LineWidth',2,'color','c','MarkerSize',20,'MarkerFaceColor','c');
    plot(range2,meanMT25fit,'.','LineWidth',2,'color','k','MarkerSize',30);
    axis([range2(1)-10,range2(end)+10,min(allPlotted25)-2,max(allPlotted25)+2])

    set(gca,'XTickMode', 'Auto','tickDir','out','xcolor','k','ycolor','k',...
        'box','off','fontname','lato','fontsize',12,'linewidth',3,...
        'ticklength',[0.02 0.02]);

    % raw data, MT-25
    retX = Cn25(n25,1);
    retX = retX(whichComparisonIndices);

    retY = Cn25(n25,2);
    retY = retY(whichComparisonIndices);

    gravX = Cn25(n25Default,1);
    gravX = gravX(whichComparisonIndices);

    gravY = Cn25(n25Default,2);
    gravY = gravY(whichComparisonIndices);

    % boxcar smoothed data, MT-25
    subplot(2,1,1)
    hold on;

    yyaxis left
    barHandles = length(range2);

    for rangeIdx = 1:1:length(range2)
        barHandles(rangeIdx) = bar(range2(rangeIdx),meanMTn25fit(rangeIdx),'BarWidth', 5,'LineWidth',0.1,'EdgeColor','w');
        set(barHandles(rangeIdx), 'FaceColor', [0.9,0.9,0.9],'EdgeColor','none');
    end

    axis([range2(1)-10,range2(end)+10,min(allPlottedn25)-2,max(allPlottedn25)+2])
    xlabel('Object Lean in Degrees','fontname','lato','FontSize', 15)
    ylabel('Firing Rate','fontname','lato','FontSize', 15)
    title('-25º Monkey Tilt, Boxcar Smoothed Data, 0º Horizon','fontname','lato','FontSize',15)

    yyaxis right
    plot(Cn25(:,1),Cn25(:,2),'LineWidth',2,'color','k');
    plot(gravX,gravY,'v','LineWidth',2,'color','m','MarkerSize',20,'MarkerFaceColor','m');
    plot(retX,retY,'^','LineWidth',2,'color','c','MarkerSize',20,'MarkerFaceColor','c');
    plot(range2,meanMTn25fit,'.','LineWidth',2,'color','k','MarkerSize',30);
    axis([range2(1)-10,range2(end)+10,min(allPlottedn25)-2,max(allPlottedn25)+2])

    set(gca,'XTickMode', 'Auto','tickDir','out','xcolor','k','ycolor','k',...
        'box','off','fontname','lato','fontsize',12,'linewidth',3,...
        'ticklength',[0.02 0.02]);

    set(gcf, 'units','normalized','outerposition',[0 0 0.5 1]);

    close all;
end

%%

function C = createCatmullRomSpline(points,context)
    curveExtent = length(points);
    C = points(1,:);

    for pointIdx = 1:1:curveExtent-3
        workingPts = points(pointIdx:pointIdx+3,:);
        C = cat(1,C,calculateSpline(workingPts,context));
    end

    C = cat(1,C,points(end,:));
end

function C = calculateSpline(workingPts,context)
    ts = [0];

    for getT = 1:1:3
        % get t1, t2, t3 beginning with t0 = 0
        tnext = findTnext(ts(getT),workingPts(getT,:),workingPts(getT+1,:));
        ts = cat(1,ts,tnext);
    end
    
    t0 = ts(1);
    t1 = ts(2);
    t2 = ts(3);
    t3 = ts(4);
    
    if contains(context,'original')
        nPoints = 26;

    elseif contains(context,'revised')
        nPoints = 51;
    end
    
    t = linspace(t1,t2,nPoints);
    
    A1 = ((t1-t)./(t1-t0))'.*workingPts(1,:)+((t-t0)./(t1-t0))'.*workingPts(2,:);
    A2 = ((t2-t)./(t2-t1))'.*workingPts(2,:)+((t-t1)./(t2-t1))'.*workingPts(3,:);
    A3 = ((t3-t)./(t3-t2))'.*workingPts(3,:)+((t-t2)./(t3-t2))'.*workingPts(4,:);
    
    B1 = ((t2-t)./(t2-t0))'.*A1+((t-t0)./(t2-t0))'.*A2;
    B2 = ((t3-t)./(t3-t1))'.*A2+((t-t1)./(t3-t1))'.*A3;
    
    C = ((t2-t)/(t2-t1))'.*B1+((t-t1)./(t2-t1))'.*B2;
end

function tnext = findTnext(tnow,pointnow,pointnext)
    alpha = 0.5;
    tnext = sqrt((pointnext(1)-pointnow(1))^2+(pointnext(2)-pointnow(2))^2)^alpha + tnow;
end
