
close all;
clear all;

load('figure4_horizonTilt.mat');

%%

colorScheme = [236,236,237;0,214,214;255,64,141]./255;
both = [155,176,44]./255;

%% HORIZON TILT SCATTER HIST FULL RETINAL

cellCorrGravitationalTemp = [];
cellCorrRetinalTemp = [];
whichSignificant = [];

for n = 1:1:length(cellPearsonCorrRetinalFullRet25n25)
    
    if significancePearsonCorrGravitational25n25(n)&&significancePearsonCorrFullRetinal25n25(n)
        whichSignificant = cat(1,whichSignificant,"rgboth");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinalFullRet25n25(n));
        
        whichSignificant = cat(1,whichSignificant,"gravitational");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinalFullRet25n25(n));
        
        whichSignificant = cat(1,whichSignificant,"retinal");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinalFullRet25n25(n));
        
    elseif significancePearsonCorrGravitational25n25(n)
        whichSignificant = cat(1,whichSignificant,"gravitational");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinalFullRet25n25(n));
        
    elseif significancePearsonCorrFullRetinal25n25(n)
        whichSignificant = cat(1,whichSignificant,"retinal");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinalFullRet25n25(n));
    end
end
    
additionalSpecification = significancePearsonCorrFullRetinal25n25|significancePearsonCorrRetinal25n25;
allGrays = ones(length(cellPearsonCorrRetinalFullRet25n25),1);

cellCorrGravitationalTempCat = cat(1,cellPearsonCorrGravitational25n25,cellPearsonCorrGravitational25n25(additionalSpecification));
cellCorrRetinalTempCat = cat(1,cellPearsonCorrRetinalFullRet25n25,cellPearsonCorrRetinalFullRet25n25(additionalSpecification));
whichSignificantCat = cat(1,allGrays,allGrays(additionalSpecification).*5);

cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTempCat,cellCorrGravitationalTemp);
cellCorrRetinalTemp = cat(1,cellCorrRetinalTempCat,cellCorrRetinalTemp);
whichSignificant = cat(1,whichSignificantCat,whichSignificant);

colorsScatterhist(1,:) = colorScheme(1,:); % all, light gray
colorsScatterhist(2,:) = [0.5,0.5,0.5]; % select ret, dark gray
colorsScatterhist(3,:) = colorScheme(3,:); % grav, pink
colorsScatterhist(4,:) = colorScheme(2,:); % full ret, teal
colorsScatterhist(5,:) = both; % both

gcf = figure;
h = scatterhist(cellCorrRetinalTemp,cellCorrGravitationalTemp,'Group',whichSignificant,'Color',colorsScatterhist,'Kernel','off');

colorsScatterhist(1,:) = colorScheme(1,:); % all, light gray
colorsScatterhist(2,:) = [0.5,0.5,0.5]; % select ret, dark gray
colorsScatterhist(3,:) = colorScheme(3,:); % grav, pink
colorsScatterhist(4,:) = colorScheme(2,:); % full ret, teal
colorsScatterhist(5,:) = both; % both

cla(h(2))
g = findgroups(whichSignificant); 
hold(h(2),'on')
arrayfun(@(i)histogram(h(2),cellCorrRetinalTemp(g==i),'BinWidth',0.05,'FaceColor',colorsScatterhist(i,:),'EdgeColor','none','FaceAlpha',1),unique(g))

set(h(2),'TickDir','out','xcolor','k','ycolor','k',...
                'box','on','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02],'xlim',[-1.1,1.1],'ylim',[0,20])

cla(h(3))
hold(h(3),'on')
arrayfun(@(i)histogram(h(3),cellCorrGravitationalTemp(g==i),'BinWidth',0.05,'FaceColor',colorsScatterhist(i,:),'EdgeColor','none','FaceAlpha',1),unique(g))

set(h(3),'TickDir','out','xcolor','k','ycolor','k',...
                'box','on','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02],'xlim',[-1.1,1.1],'ylim',[0,20])
            
set(gcf,'units','normalized','outerposition',[0 0 0.5 1.0]);
axis(h(2:3), 'on')

axis([-1.1,1.1,-1.1,1.1])
axis square;

title('Colored By Significance of Pearson Correlation, Horizon Tilt [Figure 4A]','fontname','lato','FontSize', 10)
xlabel('Retinal Pearson Correlation Coefficient','fontname','lato','FontSize', 10)
ylabel('Gravitational Pearson Correlation Coefficient','fontname','lato','FontSize', 10)
legend({'All Cells','Select Retinal','Gravitational','Full Retinal','Both'},'FontName','lato','FontSize',8,'Box','off','location','southwest')

set(gca,'TickDir','out','xcolor','k','ycolor','k',...
                'box','off','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02])

            
%% HORIZON TILT SCATTER HIST SELECT RETINAL

cellCorrGravitationalTemp = [];
cellCorrRetinalTemp = [];
whichSignificant = [];

for n = 1:1:length(cellPearsonCorrRetinal25n25)
    
    if significancePearsonCorrGravitational25n25(n)&&significancePearsonCorrRetinal25n25(n)
        whichSignificant = cat(1,whichSignificant,"rgboth");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal25n25(n));
        
        whichSignificant = cat(1,whichSignificant,"gravitational");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal25n25(n));
        
        whichSignificant = cat(1,whichSignificant,"retinal");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal25n25(n));
        
    elseif significancePearsonCorrGravitational25n25(n)
        whichSignificant = cat(1,whichSignificant,"gravitational");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal25n25(n));
        
    elseif significancePearsonCorrRetinal25n25(n)
        whichSignificant = cat(1,whichSignificant,"retinal");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational25n25(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal25n25(n));
    end
end
    
additionalSpecification = significancePearsonCorrFullRetinal25n25|significancePearsonCorrRetinal25n25;
allGrays = ones(length(cellPearsonCorrRetinal25n25),1);

cellCorrGravitationalTempCat = cat(1,cellPearsonCorrGravitational25n25,cellPearsonCorrGravitational25n25(additionalSpecification));
cellCorrRetinalTempCat = cat(1,cellPearsonCorrRetinal25n25,cellPearsonCorrRetinal25n25(additionalSpecification));
whichSignificantCat = cat(1,allGrays,allGrays(additionalSpecification).*5);

cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTempCat,cellCorrGravitationalTemp);
cellCorrRetinalTemp = cat(1,cellCorrRetinalTempCat,cellCorrRetinalTemp);
whichSignificant = cat(1,whichSignificantCat,whichSignificant);

colorsScatterhist(1,:) = colorScheme(1,:); % all, light gray
colorsScatterhist(2,:) = [0.5,0.5,0.5]; % select ret, dark gray
colorsScatterhist(3,:) = colorScheme(2,:); % full ret, teal
colorsScatterhist(4,:) = colorScheme(3,:); % grav, pink
colorsScatterhist(5,:) = both; % both

gcf = figure;
h = scatterhist(cellCorrRetinalTemp,cellCorrGravitationalTemp,'Group',whichSignificant,'Color',colorsScatterhist,'Kernel','off');

colorsScatterhist(1,:) = colorScheme(1,:); % all, light gray
colorsScatterhist(2,:) = [0.5,0.5,0.5]; % select ret, dark gray
colorsScatterhist(3,:) = colorScheme(3,:); % grav, pink
colorsScatterhist(4,:) = colorScheme(2,:); % full ret, teal
colorsScatterhist(5,:) = both; % both

cla(h(2))
g = findgroups(whichSignificant); 
hold(h(2),'on')
arrayfun(@(i)histogram(h(2),cellCorrRetinalTemp(g==i),'BinWidth',0.05,'FaceColor',colorsScatterhist(i,:),'EdgeColor','none','FaceAlpha',1),unique(g))

set(h(2),'TickDir','out','xcolor','k','ycolor','k',...
                'box','on','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02],'xlim',[-1.1,1.1],'ylim',[0,20])

cla(h(3))
hold(h(3),'on')
arrayfun(@(i)histogram(h(3),cellCorrGravitationalTemp(g==i),'BinWidth',0.05,'FaceColor',colorsScatterhist(i,:),'EdgeColor','none','FaceAlpha',1),unique(g))

set(h(3),'TickDir','out','xcolor','k','ycolor','k',...
                'box','on','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02],'xlim',[-1.1,1.1],'ylim',[0,20])
            
set(gcf,'units','normalized','outerposition',[0 0 0.5 1.0]);
axis(h(2:3), 'on')

axis([-1.1,1.1,-1.1,1.1])
axis square;

title('Colored By Significance of Pearson Correlation, Horizon Tilt [Figure 4B]','fontname','lato','FontSize', 10)
xlabel('Retinal Pearson Correlation Coefficient','fontname','lato','FontSize', 10)
ylabel('Gravitational Pearson Correlation Coefficient','fontname','lato','FontSize', 10)
legend({'All Cells','Select Retinal','Full Retinal','Gravitational','Both'},'FontName','lato','FontSize',8,'Box','off','location','southwest')

set(gca,'TickDir','out','xcolor','k','ycolor','k',...
                'box','off','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02])

            
            
            
            