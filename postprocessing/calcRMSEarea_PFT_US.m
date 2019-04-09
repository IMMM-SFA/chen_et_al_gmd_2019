function calcRMSEarea_PFT_US(year)

dir = '/pic/projects/im3/lulc_thrust_area/minchen/demeter_calibration/ensemble/';
%year = 2015;
refdata = importdata([dir '/tmp/refdata/ESACCI-LC-aggregated-0.25Deg-merged' num2str(year) '-reordered.csv']);
refdata = refdata.data;
lat = refdata(:,1);
lon = refdata(:,2);

regfilter = lat < 50 & lat > 25 & lon< -65 & lon > -125;
area = importdata('demeter_grid_area.mat');
area = area(regfilter);

RMSE = NaN(23100,7);

parpool(12)
parfor j=1:23100
    j
    fn = [dir 'outputs/outputs_' num2str(j) '/outputs/outputs_' num2str(j) '/spatial_landcover_tabular/landcover_' num2str(year) '_timestep.csv'];
    tmp = importdata(fn);
        
    for i = 1:7
        data = tmp.data(regfilter,i+1);
        ref = refdata(regfilter,i+2);

        RMSE(j,i) = sqrt(nanmean(data.*area-ref.*area).^2);
    end
end

save(['RMSEarea_PFT_US_' num2str(year) '.mat'], 'RMSE','-v7.3');
