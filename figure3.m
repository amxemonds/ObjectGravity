
close all;
clear all;

load('horizonTilt_catmullData.mat')

%% TILTING HORIZONS

for cell = 117:1:numCells

    if cell < 48
        range2 = [-50,-37.5,-25,-12.5,0,12.5,25,37.5,50];   % (50 to -50)
        rangeObjOri = 1:9;
        matchProbe = 5;
        catromRef = 62.5;
        context = 'original';

    elseif cell <= 118
        range2 = [-100,-75,-50,-25,0,25,50,75,100];         % (100 to -100)
        rangeObjOri = 2:10;
        matchProbe = 7;
        catromRef = 125;
        context = 'revised';
    else
        range2 = [-50,-37.5,-25,-12.5,0,12.5,25,37.5,50];   % (50 to -50)
        rangeObjOri = 1:9;
        matchProbe = 5;
        catromRef = 62.5;
        context = 'original';
    end

    % these are labeled -50 to 50, but are actually 50 to -50
    % do flipud to correct here

    % have one distribution for HT25
    meanHT25raw = recordAccumulationCellsP25(:,1,cell);
    meanHT25raw(rangeObjOri) = flipud(meanHT25raw(rangeObjOri));

    % have one distribution for HT-25
    meanHTn25raw = recordAccumulationCellsN25(:,1,cell);
    meanHTn25raw(rangeObjOri) = flipud(meanHTn25raw(rangeObjOri));

    if (sum(meanHT25raw)==0) || (sum(meanHTn25raw)==0)
        continue;
    end

    % find the values at the smooth fit for HT25
    meanHT25rawPadded = padarray(meanHT25raw(rangeObjOri),[1,0],'replicate');
    meanHT25fit = conv(meanHT25rawPadded',[1/3,1/3,1/3],'valid')';

    % find the values at the smooth fit for HT-25
    meanHTn25rawPadded = padarray(meanHTn25raw(rangeObjOri),[1,0],'replicate');
    meanHTn25fit = conv(meanHTn25rawPadded',[1/3,1/3,1/3],'valid')';

    yQuarticN = meanHTn25fit;
    yQuarticP = meanHT25fit;

    % catmull-rom splines for shifted analysis
    rangePadded = [-catromRef,range2,catromRef]';

    meanHT25fitPadded = padarray(meanHT25fit,[1,0],'replicate');
    C25 = createCatmullRomSpline([rangePadded,meanHT25fitPadded],context);
    C25 = C25(2:end-1,:);

    for nn = 2:length(range2)-1
        idx = find(round(C25(:,1),1)==range2(nn));

        if idx
            C25(idx(1),:) = []; 
        end
    end

    C25 = unique(C25,'rows');

    meanHTn25fitPadded = padarray(meanHTn25fit,[1,0],'replicate');
    Cn25 = createCatmullRomSpline([rangePadded,meanHTn25fitPadded],context);
    Cn25 = Cn25(2:end-1,:);

    for nn = 2:length(range2)-1
        idx = find(round(Cn25(:,1),1)==range2(nn));

        if idx
            Cn25(idx(1),:) = []; 
        end
    end

    Cn25 = unique(Cn25,'rows');

    gcf = figure;
    allPlotted25 = [meanHT25raw(rangeObjOri);meanHT25fit];
    allPlottedn25 = [meanHTn25raw(rangeObjOri);meanHTn25fit];

    % raw data, HT25
    if contains(context,'original')
        gravX = range2(5:9);
        gravY = yQuarticP(5:9);

        retX = range2(5:9);
        retY = yQuarticP(5:9);

    elseif contains(context,'revised')
        gravX = range2(3:9);
        gravY = yQuarticP(3:9);

        retX = range2(3:9);
        retY = yQuarticP(3:9);
    end

    % boxcar smoothed data, MT25
    subplot(2,1,2)
    hold on;

    yyaxis left
    barHandles = length(range2);

    for rangeIdx = 1:1:length(range2)
        barHandles(rangeIdx) = bar(range2(rangeIdx),meanHT25fit(rangeIdx),'BarWidth', 5,'LineWidth',0.1,'EdgeColor','w');
        set(barHandles(rangeIdx), 'FaceColor', [0.9,0.9,0.9],'EdgeColor','none');
    end

    axis([range2(1)-10,range2(end)+10,min(allPlotted25)-2,max(allPlotted25)+2])
    xlabel('Object Lean in Degrees','fontname','lato','FontSize', 15)
    ylabel('Firing Rate','fontname','lato','FontSize', 15)
    title('+25º Horizon Tilt, Boxcar Smoothed Data','fontname','lato','FontSize',15)

    yyaxis right
    plot(C25(:,1),C25(:,2),'LineWidth',2,'color','k');
    plot(retX,retY,'^','LineWidth',2,'color','c','MarkerSize',20,'MarkerFaceColor','c');
    plot(gravX,gravY,'v','LineWidth',2,'color','m','MarkerSize',20,'MarkerFaceColor','m');
    plot(range2,meanHT25fit,'.','LineWidth',2,'color','k','MarkerSize',30);
    axis([range2(1)-10,range2(end)+10,min(allPlotted25)-2,max(allPlotted25)+2])

    set(gca,'XTickMode', 'Auto','tickDir','out','xcolor','k','ycolor','k',...
        'box','off','fontname','lato','fontsize',12,'linewidth',3,...
        'ticklength',[0.02 0.02]);

    % raw data, MT-25
    if contains(context,'original')
        gravX = range2(1:5);
        gravY = yQuarticN(1:5);

        retX = range2(5:9);
        retY = yQuarticN(5:9);

    elseif contains(context,'revised')
        gravX = range2(1:7);
        gravY = yQuarticN(1:7);

        retX = range2(3:9);
        retY = yQuarticN(3:9);
    end

    % boxcar smoothed data, MT-25
    subplot(2,1,1)
    hold on;

    yyaxis left
    barHandles = length(range2);

    for rangeIdx = 1:1:length(range2)
        barHandles(rangeIdx) = bar(range2(rangeIdx),meanHTn25fit(rangeIdx),'BarWidth', 5,'LineWidth',0.1,'EdgeColor','w');
        set(barHandles(rangeIdx), 'FaceColor', [0.9,0.9,0.9],'EdgeColor','none');
    end

    axis([range2(1)-10,range2(end)+10,min(allPlottedn25)-2,max(allPlottedn25)+2])
    xlabel('Object Lean in Degrees','fontname','lato','FontSize', 15)
    ylabel('Firing Rate','fontname','lato','FontSize', 15)
    title('-25º Horizon Tilt, Boxcar Smoothed Data','fontname','lato','FontSize',15)

    yyaxis right
    plot(Cn25(:,1),Cn25(:,2),'LineWidth',2,'color','k');
    plot(retX,retY,'^','LineWidth',2,'color','c','MarkerSize',20,'MarkerFaceColor','c');
    plot(gravX,gravY,'v','LineWidth',2,'color','m','MarkerSize',20,'MarkerFaceColor','m');
    plot(range2,meanHTn25fit,'.','LineWidth',2,'color','k','MarkerSize',30);
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

