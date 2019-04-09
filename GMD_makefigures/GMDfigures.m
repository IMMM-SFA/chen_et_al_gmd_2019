%% add mapping tools

addpath('m_map');

%% Read in the errors of the ensemble downscaling

years = [1992,2000,2005,2010,2015];
vtypes = {'(a) Forest','(b) Shrub','(c) Grass','(d) Crop','(e) Urban','(f) Snow','(g) Sparse','(h) All FLTs'};

E = NaN(6,8,23100);
for i = 1:length(years)
    year = years(i);
    fn = ['RMSEarea_PFT_year_' num2str(year) '.mat'];
    tmpE = importdata(fn);
    E(i,1:7,:) = tmpE;
    
    % all FLTs mean
    E(i,8,:) = mean(tmpE);
end
% mean of all years
E(6,:,:) = squeeze(mean(E(1:5,:,:)));

clear tmpE;

%% Figure 2

AE = squeeze(E(6,8,:));
[b,I] = sort(AE);
allpara = importdata('allparam.csv');

figure;

subplot(3,1,1)
histogram(b,'Normalization','probability')
xlim([8,17]);
ylim([0,0.1]);
hold on;
l1 = plot(prctile(b,5)*ones(1,11),0:0.01:0.1,'k--','linewidth',2); 
l2 = plot(prctile(b,10)*ones(1,11),0:0.01:0.1,'k:','linewidth',2); 
xlabel('E (Unit: km^2)');
ylabel('Probability');
legend([l1,l2], {'5% interval','10% interval'});
legend boxoff;

set(gca,'FontSize',20)
textLoc('(a)','NorthWest','FontSize',20,'Color','k');

subplot(3,1,2)
worsti = I(b>prctile(AE,95));
besti = I(b<prctile(AE,5));
g = allpara.data(:,2:end);
g(:,6) = g(:,6)/100;
gworst = g(worsti,:);
gbest = g(besti,:);

for i = 1:length(gbest)
    plot(gbest(i,:),'-','color',[0.5,0.5,0.5],'linewidth',1);
    hold on;
end
h=violinplot(gbest,{'w_N','w_S','w_K','r','\tau', 'D'},'Width',0.3,'ShowData',false,'ShowMean',true);
plot(g(I(1),:),'r-','linewidth',3);
textLoc('(b)','NorthWest','FontSize',20,'Color','k');
set(gca,'FontSize',20)
xlabel('Demeter Parameters');
ylabel('Parameter values (X100 for D)');

subplot(3,1,3)
besti = I(b<prctile(AE,10));
g = allpara.data(:,2:end);
g(:,6) = g(:,6)/100;
gbest = g(besti,:);

for i = 1:length(gbest)
    plot(gbest(i,:),'-','color',[0.5,0.5,0.5],'linewidth',1);
    hold on;
end
h=violinplot(gbest,{'w_N','w_S','w_K','r','\tau', 'D'},'Width',0.3,'ShowData',false,'ShowMean',true);
plot(g(I(1),:),'r-','linewidth',3);
textLoc('(c)','NorthWest','FontSize',20,'Color','k');
set(gca,'FontSize',20)
xlabel('Demeter Parameters');

%% Figure 3
allpara = importdata('allparam.csv');
g = allpara.data(:,2:end);

paranames = {'w_N','w_S','w_K','r','\tau', 'D'};
abcd = 'abcdefg';
figure;
for i = 1:6
    subplot(2,3,i)
    
    gi = g(:,i);
    Ei = squeeze(E(end,end,:));
        
    boxplot(Ei,gi,'Colors',[0,0.447,0.741],...
        'Widths',0.5,'BoxStyle','outline',...
        'MedianStyle','target','Symbol','k.',...
        'PlotStyle','traditional');
    hold on;
    
    b = polyfit(gi,Ei,1)
    
    ugi = unique(gi);
    x = 1:length(ugi);
    plot(x,ugi(x)*b(1)+b(2),'k-','linewidth',3);
    
    [rho,pval] = corr(gi,Ei);
    set(gca,'FontSize',25);
    
    ylim([8 ,18]);
    xlabel(paranames{i});
    ylabel('E (km^2)')
    % txt = ['R^2=' num2str(rho^2,'%2.3f') '; p=' num2str(pval,'%2.2f')];
    txt = ['R^2=' num2str(rho^2,'%2.3f') '; p<0.01'];
    textLoc(txt,'NorthEast','FontSize',25,'Color','k');
    textLoc(['(' abcd(i) ')'],'NorthWest','FontSize',25,'Color','k');
    
end

%% Figure 4
% the Sobol indices were calculated with in R with package "sensitivity"

paranames = {'w_N','w_S','w_K','r','\tau', 'D'};

SI1 = [ 1.672004e-03  0.02590065;
    4.204563e-03  0.02387416;
    1.236419e-05  0.01440629;
    6.785759e-01  0.76760898;
    2.063360e-01  0.30728261;
    4.352264e-03  0.02095197];

SI2 = [ 0.0156055933  0.06643164;
    0.0025121890  0.04547865;
    0.0000641747  0.03716257;
    0.3094969567  0.43560930;
    0.5095542926  0.64770635;
    0.0141776415  0.05678577];

SI3 = [ 1.154994e-04  0.01510649;
    1.967032e-05  0.01459152;
    1.090052e-06  0.01451717;
    6.030961e-01  0.71083211;
    2.568158e-01  0.36310021;
    2.402956e-02  0.05341455];

SI4 = [ 3.648047e-03  0.04383648;
    5.899818e-03  0.03863573;
    2.921966e-05  0.02096156;
    6.002713e-01  0.71022262;
    2.529118e-01  0.38581682;
    5.007081e-04  0.02609639];

SI5 = [1.118805e-06  0.01150447;
    5.076097e-06  0.01149877;
    4.348108e-07  0.01147631;
    4.542960e-01  0.59044540;
    3.802181e-01  0.51792080;
    2.662782e-02  0.04296018];

SI6 = [ 3.676810e-04  0.07480222;
    3.796945e-04  0.07483272;
    1.427426e-06  0.07356695;
    8.079416e-02  0.34469808;
    6.317096e-01  0.89883330;
    8.321158e-03  0.11395870];

SI7 = [ 1.426716e-03  0.02729870;
    1.420101e-03  0.02628360;
    9.982956e-06  0.01944049;
    5.643642e-01  0.69292797;
    2.868422e-01  0.42320202;
    8.814055e-03  0.03859704];

SIall = [ 1.327412e-03  0.01996015;
    1.531851e-03  0.01870720;
    9.370905e-06  0.01428120;
    5.886612e-01  0.69275060;
    2.860539e-01  0.39886341;
    8.751193e-03  0.02694141];

vtypes = {'(a) Forest','(b) Shrub','(c) Grass','(d) Crop','(e) Urban','(f) Snow','(g) Sparse','(h) All FLTs'};

figure;

for i=1:7
    subplot(2,4,i)
    if i <= 7
        eval(['SI = SI' num2str(i) ';']);
    else
        SI = SIall;
    end
    bar(1:6,SI);
    colormap('Gray')
    set(gca, 'XTickLabel', paranames,'FontSize',30)
    xlabel('Parameters');
    ylabel('Sobol Indices');
    xlim([0,7]);
    ylim([0,1]);
    if i==1
        legend('First order index','Total order index');
    end
    title(vtypes{i})
end

%% Figure 5

figure;
for i = 1:8
    tmpE = squeeze(E(:,i,:))';
    
    for j = 1:6
        if std(tmpE(:,j)) < mean(tmpE(:,j))/10
        res = rand(size(tmpE(:,j)))*mean(tmpE(:,j))/10;
        tmpE(:,j) = tmpE(:,j)+res;
        end
    end
    subplot(2,4,i)
    h=violinplot(tmpE,{'1992','2000','2005','2010','2015','Mean'},'Width',0.3,'ShowData',false,'ShowMean',true);
    xtickangle(45);
    axis square;
    set(gca,'FontSize',35);
    ylim([0,45]);
    xlim([0,7]);
    ylabel('E (km^2)')
    box on;
    title(vtypes{i})
end


%% Figure 6

refdata = importdata(['ESACCI-LC-aggregated-0.25Deg-merged2015-reordered.csv']);
refdata = refdata.data;
lat = refdata(:,1);
lon = refdata(:,2);
bestresult = importdata('downscaled_4360.csv');
bestresult = bestresult.data;

vtypes = {'(a) Forest','(b) Shrub','(c) Grass','(d) Crop','(e) Urban','(f) Snow','(g) Sparse','(h) All FLTs'};

figure;
for i = 1:7
    subplot(2,4,i)
    x = refdata(:,2+i)*100;
    y = bestresult(:,1+i)*100;
    r2 = corr(x,y)^2;
    plot(x,y,'k.');
    hold on;
    
    l = refline(1,0);
    set(l,'linewidth',2.5,'color','b');
    
    [p,S] = polyfit(x,y,1);
    xfit = 0:100;
    [Y,DELTA] = polyconf(p,xfit,S,'alpha',0.05);
    plot(xfit,Y+DELTA,'r--','linewidth',2.5);
    plot(xfit,Y-DELTA,'r--','linewidth',2.5);
    
    
    xlim([0,100]);
    ylim([0,100]);
    str = ['R^2=' num2str(r2,'%0.2f')];
    textLoc(str,'NorthWest','FontSize',30,'FontWeight','bold','Color','k');
    title(vtypes{i});
    set(gca,'FontSize',30);
    xlabel('Fraction of grid (%)');
    ylabel('Fraction of grid (%)');
end

%% Figure 7

% change the number for other FLTs as the figures in the supplementary
% materials Figure S2-S6
iFLT = 7;

vtypes = {'Forest','Shrub','Grass','Crop','Urban','Snow','Sparse','All FLTs'};

reffrac = refdata(:,iFLT+2);
bestfrac = bestresult(:,iFLT+1);

refmap = makemap(lat,lon,reffrac,0.25)*100;
bestmap = makemap(lat,lon,bestfrac,0.25)*100;
mapdif = bestmap-refmap;

M=m_shaperead('aezshapefiles/aez_orig_lds_5arcmin');
[Mlon,Mlat] = meshgrid(-179.875:0.25:179.875,89.875:-0.25:-89.875);

colorn = 'RdBu11';
cmap = othercolor(colorn,1000);
cmap = flipud(cmap(1:500,:));

figure;
subplot(3,1,1)
m_proj('Robinson','lon',[-180 180],'lat',[-90 90]);
hold on;
m_pcolor(Mlon,Mlat,refmap);
shading interp;
for k=1:length(M.ncst)
    m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'linestyle',':','linewidth',0.1,'color',[0.6 0.6 0.6]);
end
m_coast('line','color','k');
m_grid('box','on','tickdir','in');
colormap(cmap);
caxis([0,100]);
set(gca,'FontSize',18);
hcb=colorbar;
hcb.Label.String = [vtypes{iFLT} ' fraction - reference (%)'];
textLoc('(a)','NorthWest','FontSize',18,'Color','k');

subplot(3,1,2)
m_proj('Robinson','lon',[-180 180],'lat',[-90 90]);
hold on;
m_pcolor(Mlon,Mlat,bestmap);
shading interp;
for k=1:length(M.ncst)
    m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'linestyle',':','linewidth',0.1,'color',[0.6 0.6 0.6]);
end
m_coast('line','color','k');
m_grid('box','on','tickdir','in');
colormap(cmap);
caxis([0,0100]);
set(gca,'FontSize',18);
hcb=colorbar;
hcb.Label.String = [vtypes{iFLT} 'fraction - downscaled (%)'];
textLoc('(b)','NorthWest','FontSize',18,'Color','k');

ax3=subplot(3,1,3)
colorn = 'RdBu11';
cmap = othercolor(colorn,1000);

m_proj('Robinson','lon',[-180 180],'lat',[-90 90]);
hold on;
m_pcolor(Mlon,Mlat,mapdif);
shading interp;
for k=1:length(M.ncst)
    m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'linestyle',':','linewidth',0.1,'color',[0.6 0.6 0.6]);
end
m_coast('line','color','k');
m_grid('box','on','tickdir','in');
colormap(ax3,cmap);
caxis([-10,10]);
set(gca,'FontSize',18);
hcb=colorbar;
hcb.Label.String = 'Difference (downscaled-reference) (%)';
textLoc('(c)','NorthWest','FontSize',18,'Color','k');

%% Figure 8
area = importdata('area.mat');
area = repmat(area,1,12);
vtypes = {'Forest','Shrub','Grass','Crop','Urban','Snow','Sparse','All FLTs'};
years = [2005,2010,2050,2100];

nr = 2;
nc = 4;
outerspacer = 0.1;
outerspacec = 0.1;
marginc = 0.03;
marginr = 0.05;
width = (1-marginc*(nc-1)-outerspacec*2)/nc;
height = (1-marginr*(nr-1)-outerspacer*2)/nr;
pos=zeros(nr*nc,4);
for i=1:nr*nc
    [col,row] = ind2sub([nc nr],i);
    left = (col-1)*width+(col-1)*marginc+outerspacec;
    bott = 1-outerspacer-(row*height+(row-1)*marginr);
    pos(i,:) = [left bott width height];
end
abcd='abcdefghijklmnopqrstuvwxyz';

figure;
for iFLT = 1:7
    for i = 1:length(years)        
        ystr = num2str(years(i));
        fn = ['LC_' ystr '_std_5p.mat'];
        eval(['LCstd5p' ystr '= importdata(fn);']);
        eval(['LCstd5p' ystr '= area.*LCstd5p' ystr ';']);
        fn = ['LC_' ystr '_std_10p.mat'];
        eval(['LCstd10p' ystr '= importdata(fn);']);
        eval(['LCstd10p' ystr '= area.*LCstd10p' ystr ';']);
    end
    
    hf(iFLT) = subplot(2,4,iFLT);
    
    ind = iFLT + 5;
    
    LCstd = [LCstd5p2005(:,ind) LCstd5p2010(:,ind) LCstd5p2050(:,ind) LCstd5p2100(:,ind)];
    shadedErrorBar(years,nanmean(LCstd),nanstd(LCstd),{'-','LineWidth',3,'Color','k'}, 0.5);
    
    hold on;
    
    LCstd = [LCstd10p2005(:,ind) LCstd10p2010(:,ind) LCstd10p2050(:,ind) LCstd10p2100(:,ind)];
    shadedErrorBar(years,nanmean(LCstd),nanstd(LCstd),{'-','LineWidth',3,'Color','r'}, 0.5);
    
    ylim([0,100]);
    set(gca,'FontSize',15);
    xlabel('Year','FontSize',25);
    
    if iFLT == 1 | iFLT == 5
        ylabel('\sigma (km^2)','FontSize',25);
    end
    xticks([2005 2010 2020:10:2100]);
    xtickangle(45);
    xlim([2005,2100]);
    title(['(' abcd(iFLT) ') ' vtypes{iFLT}],'FontSize',25);  
    axis square;
end


for i=1:7
    set(hf(i),'Position',pos(i,:));
end


%% Figure S7
vtypes = {'Forest','Shrub','Grass','Crop','Urban','Snow','Sparse','All FLTs'};

colorn = 'RdBu11';
cmap = othercolor(colorn,1000);
cmap = flipud(cmap(1:500,:));

[Mlon,Mlat] = meshgrid(-179.875:0.25:179.875,89.875:-0.25:-89.875);

nr = 7;
nc = 4;
outerspacer = 0.1;
outerspacec = 0.1;
marginc = 0.00;
marginr = 0.02;
width = (1-marginc*(nc-1)-outerspacec*2)/nc;
height = (1-marginr*(nr-1)-outerspacer*2)/nr;
pos=zeros(nr*nc,4);
for i=1:nr*nc
    [col,row] = ind2sub([nc nr],i);
    left = (col-1)*width+(col-1)*marginc+outerspacec;
    bott = 1-outerspacer-(row*height+(row-1)*marginr);
    pos(i,:) = [left bott width height];
end
abcd='abcdefghijklmnopqrstuvwxyz';

figure;

years = [2005 2010 2050 2100];
k = 0;
for iFLT = 1:7
    ind = iFLT + 5;
    for i = 1:length(years)
        
        k = k+1;
        
        ystr = num2str(years(i));
        fn = ['LC_' ystr '_mean_5p.mat'];
        eval(['LCavg' ystr '= importdata(fn);']);
        fn = ['LC_' ystr '_std_5p.mat'];
        eval(['LCstd' ystr '= importdata(fn);']);
        
        eval(['LCavg = LCavg' ystr ';']);
        eval(['LCstd = LCstd' ystr '.*area;']);
        
        map1 = makemap(lat,lon,LCavg(:,ind),0.25)*100;
        map2 = makemap(lat,lon,LCstd(:,ind),0.25);
        
        hf(k) = subplot(nr,nc,k);
        surf(Mlon,Mlat,map2,map1,'EdgeColor','none');
        view(45,45);
        colormap(cmap);
        caxis([0,100]);
        xlim([-180,180]);
        ylim([-90,90]);
        zlim([0,400]);
        xlabel('longitude','Rotation',-38);
        ylabel('latitude','Rotation',38);
        zlabel('\sigma (km^2)');
        
        
        set(gca,'FontSize',12);
        if iFLT==1
            title(ystr,'FontSize',20);
        end
        if i==1
            textLoc(vtypes{iFLT},{'westoutside',1/2},'FontSize',20,'FontWeight','bold');
        end
        if k==12
            hcb = colorbar;
            ylabel(hcb, 'fraction of grid area (%)','FontSize',20)
        end
        axis square;
    end
end

for i=1:28
    set(hf(i),'Position',pos(i,:));
end

%% Figure S8
vtypes = {'Forest','Shrub','Grass','Crop','Urban','Snow','Sparse','All FLTs'};

colorn = 'RdBu11';
cmap = othercolor(colorn,1000);
cmap = flipud(cmap(1:500,:));

[Mlon,Mlat] = meshgrid(-179.875:0.25:179.875,89.875:-0.25:-89.875);

nr = 7;
nc = 4;
outerspacer = 0.1;
outerspacec = 0.1;
marginc = 0.00;
marginr = 0.02;
width = (1-marginc*(nc-1)-outerspacec*2)/nc;
height = (1-marginr*(nr-1)-outerspacer*2)/nr;
pos=zeros(nr*nc,4);
for i=1:nr*nc
    [col,row] = ind2sub([nc nr],i);
    left = (col-1)*width+(col-1)*marginc+outerspacec;
    bott = 1-outerspacer-(row*height+(row-1)*marginr);
    pos(i,:) = [left bott width height];
end
abcd='abcdefghijklmnopqrstuvwxyz';

figure;

years = [2005 2010 2050 2100];
k = 0;
for iFLT = 1:7
    ind = iFLT + 5;
    for i = 1:length(years)
        
        k = k+1;
        
        ystr = num2str(years(i));
        fn = ['LC_' ystr '_mean_10p.mat'];
        eval(['LCavg' ystr '= importdata(fn);']);
        fn = ['LC_' ystr '_std_10p.mat'];
        eval(['LCstd' ystr '= importdata(fn);']);
        
        eval(['LCavg = LCavg' ystr ';']);
        eval(['LCstd = LCstd' ystr '.*area;']);
        
        map1 = makemap(lat,lon,LCavg(:,ind),0.25)*100;
        map2 = makemap(lat,lon,LCstd(:,ind),0.25);
        
        hf(k) = subplot(nr,nc,k);
        surf(Mlon,Mlat,map2,map1,'EdgeColor','none');
        view(45,45);
        colormap(cmap);
        caxis([0,100]);
        xlim([-180,180]);
        ylim([-90,90]);
        zlim([0,400]);
        xlabel('longitude','Rotation',-38);
        ylabel('latitude','Rotation',38);
        zlabel('\sigma (km^2)');
        
        
        set(gca,'FontSize',12);
        if iFLT==1
            title(ystr,'FontSize',20);
        end
        if i==1
            textLoc(vtypes{iFLT},{'westoutside',1/2},'FontSize',20,'FontWeight','bold');
        end
        if k==12
            hcb = colorbar;
            ylabel(hcb, 'fraction of grid area (%)','FontSize',20)
        end
        axis square;
    end
end

for i=1:28
    set(hf(i),'Position',pos(i,:));
end

%% Figure S9
years = [1992 2000 2005 2010 2015];

E_IN = NaN(6,23100,7);
E_US = NaN(6,23100,7);
for i = 3:5
    year = years(i);
    ystr = num2str(year);
    fn = ['RMSEarea_PFT_IN_' ystr '.mat'];
    E_IN(i,:,:) = importdata(fn);
    fn = ['RMSEarea_PFT_US_' ystr '.mat'];
    E_US(i,:,:) = importdata(fn);
end
E_IN_mean = squeeze(nanmean(nanmean(E_IN),3));
E_US_mean = squeeze(nanmean(nanmean(E_US),3));


allpara = importdata('allparam.csv');
g = allpara.data(:,2:end);
g(:,6) = g(:,6)/100;

figure;
for i=1:2
    if i==1
        E = E_IN_mean;
    else
        E = E_US_mean;
    end
    [b,I] = sort(E);
    besti = I(b<prctile(E,5));
    
    gbest = g(besti,:);
  
    subplot(2,1,i)   
    h=violinplot(gbest,{'w_N','w_S','w_K','r','\tau', 'D'},'Width',0.3,'ShowData',false,'ShowMean',true);
    plot(g(I(1),:),'r-','linewidth',3);
    
    if i==1
        textLoc('(a) India','NorthWest','FontSize',20,'Color','k');
    else
        textLoc('(b) USA','NorthWest','FontSize',20,'Color','k');
    end
    set(gca,'FontSize',20)
    xlabel('Demeter Parameters');
    ylabel('Parameter values (X100 for D)');
end