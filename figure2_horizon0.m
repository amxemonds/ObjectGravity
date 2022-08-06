
close all;
clear all;

load('figure2_horizon0.mat')

%%

colorScheme = [236,236,237;0,214,214;255,64,141]./255;
both = [155,176,44]./255;

significancePearsonCorrGravitational8 = cellPearsonCorrGravitational8 > percentileSignificantGravitationalPearsonCorr8(:,2);
significancePearsonCorrFullGravitational8 = cellPearsonCorrGravitationalFullGrav8 > percentileSignificantFullGravitationalPearsonCorr8(:,2);
significancePearsonCorrRetinal8 = cellPearsonCorrRetinal8 > percentileSignificantRetinalPearsonCorr8(:,2);

%% 0º HORIZON SCATTER HIST FULL GRAVITATIONAL

cellCorrGravitationalTemp = [];
cellCorrRetinalTemp = [];
whichSignificant = [];

for n = 1:1:length(cellPearsonCorrGravitationalFullGrav8)
    
    if significancePearsonCorrFullGravitational8(n)&&significancePearsonCorrRetinal8(n)
        whichSignificant = cat(1,whichSignificant,"rgboth");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitationalFullGrav8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
        
        whichSignificant = cat(1,whichSignificant,"gravitational");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitationalFullGrav8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
        
        whichSignificant = cat(1,whichSignificant,"retinal");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitationalFullGrav8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
        
    elseif significancePearsonCorrFullGravitational8(n)
        whichSignificant = cat(1,whichSignificant,"gravitational");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitationalFullGrav8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
        
    elseif significancePearsonCorrRetinal8(n)
        whichSignificant = cat(1,whichSignificant,"retinal");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitationalFullGrav8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
    end
end
    
additionalSpecification = significancePearsonCorrFullGravitational8|significancePearsonCorrGravitational8;
allGrays = ones(length(cellPearsonCorrGravitationalFullGrav8),1);

cellCorrGravitationalTempCat = cat(1,cellPearsonCorrGravitationalFullGrav8,cellPearsonCorrGravitationalFullGrav8(additionalSpecification));
cellCorrRetinalTempCat = cat(1,cellPearsonCorrRetinal8,cellPearsonCorrRetinal8(additionalSpecification));
whichSignificantCat = cat(1,allGrays,allGrays(additionalSpecification).*5);

cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTempCat,cellCorrGravitationalTemp);
cellCorrRetinalTemp = cat(1,cellCorrRetinalTempCat,cellCorrRetinalTemp);
whichSignificant = cat(1,whichSignificantCat,whichSignificant);

%%%%% CAN'T SEE THE BOTH
colorsScatterhist(1,:) = colorScheme(1,:); % all, light gray
colorsScatterhist(2,:) = [0.5,0.5,0.5]; % select ret, dark gray
colorsScatterhist(3,:) = colorScheme(3,:); % grav, pink
colorsScatterhist(4,:) = both; % both
colorsScatterhist(5,:) = colorScheme(2,:); % full ret, teal

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
                'ticklength',[0.02 0.02],'xlim',[-1.1,1.1],'ylim',[0,10])

cla(h(3))
hold(h(3),'on')
arrayfun(@(i)histogram(h(3),cellCorrGravitationalTemp(g==i),'BinWidth',0.05,'FaceColor',colorsScatterhist(i,:),'EdgeColor','none','FaceAlpha',1),unique(g))

set(h(3),'TickDir','out','xcolor','k','ycolor','k',...
                'box','on','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02],'xlim',[-1.1,1.1],'ylim',[0,10])
            
set(gcf,'units','normalized','outerposition',[0 0 0.5 1.0]);
axis(h(2:3), 'on')

axis([-1.1,1.1,-1.1,1.1])
axis square;

title('Colored By Significance of Pearson Correlation, 0º Horizon [Figure 2A]','fontname','lato','FontSize', 10)
xlabel('Retinal Pearson Correlation Coefficient','fontname','lato','FontSize', 10)
ylabel('Gravitational Pearson Correlation Coefficient','fontname','lato','FontSize', 10)
legend({'All Cells','Select Gravitational','Full Gravitational','Both','Retinal'},'FontName','lato','FontSize',8,'Box','off','location','southwest')

set(gca,'TickDir','out','xcolor','k','ycolor','k',...
                'box','off','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02])

            
%% 0º HORIZON SCATTER HIST SELECT GRAVITATIONAL

cellCorrGravitationalTemp = [];
cellCorrRetinalTemp = [];
whichSignificant = [];

for n = 1:1:length(cellPearsonCorrGravitational8)
    
    if significancePearsonCorrGravitational8(n)&&significancePearsonCorrRetinal8(n)
        whichSignificant = cat(1,whichSignificant,"rgboth");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
        
        whichSignificant = cat(1,whichSignificant,"gravitational");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
        
        whichSignificant = cat(1,whichSignificant,"retinal");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
        
    elseif significancePearsonCorrGravitational8(n)
        whichSignificant = cat(1,whichSignificant,"gravitational");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
        
    elseif significancePearsonCorrRetinal8(n)
        whichSignificant = cat(1,whichSignificant,"retinal");
        cellCorrGravitationalTemp = cat(1,cellCorrGravitationalTemp,cellPearsonCorrGravitational8(n));
        cellCorrRetinalTemp = cat(1,cellCorrRetinalTemp,cellPearsonCorrRetinal8(n));
    end
end
    
additionalSpecification = significancePearsonCorrFullGravitational8|significancePearsonCorrGravitational8;
allGrays = ones(length(cellPearsonCorrGravitational8),1);

cellCorrGravitationalTempCat = cat(1,cellPearsonCorrGravitational8,cellPearsonCorrGravitational8(additionalSpecification));
cellCorrRetinalTempCat = cat(1,cellPearsonCorrRetinal8,cellPearsonCorrRetinal8(additionalSpecification));
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
                'ticklength',[0.02 0.02],'xlim',[-1.1,1.1],'ylim',[0,10])

cla(h(3))
hold(h(3),'on')
arrayfun(@(i)histogram(h(3),cellCorrGravitationalTemp(g==i),'BinWidth',0.05,'FaceColor',colorsScatterhist(i,:),'EdgeColor','none','FaceAlpha',1),unique(g))

set(h(3),'TickDir','out','xcolor','k','ycolor','k',...
                'box','on','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02],'xlim',[-1.1,1.1],'ylim',[0,10])
            
set(gcf,'units','normalized','outerposition',[0 0 0.5 1.0]);
axis(h(2:3), 'on')

axis([-1.1,1.1,-1.1,1.1])
axis square;

title('Colored By Significance of Pearson Correlation, 0º Horizon [Figure 2C]','fontname','lato','FontSize', 10)
xlabel('Retinal Pearson Correlation Coefficient','fontname','lato','FontSize', 10)
ylabel('Gravitational Pearson Correlation Coefficient','fontname','lato','FontSize', 10)
legend({'All Cells','Select Gravitational','Full Gravitational','Retinal','Both'},'FontName','lato','FontSize',8,'Box','off','location','southwest')

set(gca,'TickDir','out','xcolor','k','ycolor','k',...
                'box','off','fontname','lato','fontsize',10,'linewidth',3,...
                'ticklength',[0.02 0.02])


            
            
            
