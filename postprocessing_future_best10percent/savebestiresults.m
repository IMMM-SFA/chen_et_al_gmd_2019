function savebestiresults(iy)

dir = '/pic/projects/im3/lulc_thrust_area/minchen/demeter_calibration/ensemble/';
years = [2005 2010 2050 2100];

besti = importdata('besti_10percent.txt');

for i = 1:length(besti)
    j = besti(i);
    fn = [dir 'outputs_future_best10percent/outputs_' num2str(j) '/outputs_future_best10percent/outputs_' num2str(j) '/spatial_landcover_tabular/landcover_' num2str(years(iy)) '_timestep.csv']  
    tmp = importdata(fn);
    data = tmp.data;

    if i==1
       alldata = NaN(size(data,1),size(data,2),length(besti));
    end
    alldata(:,:,i) = data;
end

    datamean = nanmean(alldata,3);
    datastd = nanstd(alldata,[],3);
    
    outfn = ['LC_' num2str(years(iy)) '_mean_10p.mat']
    save(outfn,'datamean','-v7.3');
    outfn = ['LC_' num2str(years(iy)) '_std_10p.mat']
    save(outfn,'datastd','-v7.3');

    datamean = nanmean(alldata(:,:,1:1155),3);
    datastd = nanstd(alldata(:,:,1:1155),[],3);

    outfn = ['LC_' num2str(years(iy)) '_mean_5p.mat']
    save(outfn,'datamean','-v7.3');
    outfn = ['LC_' num2str(years(iy)) '_std_5p.mat']
    save(outfn,'datastd','-v7.3');
