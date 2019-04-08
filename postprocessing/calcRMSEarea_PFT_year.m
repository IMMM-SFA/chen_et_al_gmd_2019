function calcRMSEarea_PFT_year(i)

year = 1991+i

dir = '/pic/projects/im3/lulc_thrust_area/minchen/demeter_calibration/ensemble/';
refdata = importdata(['ESACCI-LC-aggregated-0.25Deg-merged' num2str(year) '-reordered.csv']);
refdata = refdata.data(:,3:9);

area = importdata('demeter_grid_area.mat');
areamat = repmat(area,1,7);

RMSE = NaN(7,23100);

parfor j=1:23100
    j
    fn = [dir 'outputs/outputs_' num2str(j) '/outputs/outputs_' num2str(j) '/spatial_landcover_tabular/landcover_' num2str(year) '_timestep.csv'];
    tmp = importdata(fn);
    data = tmp.data(:,2:8);
    RMSE(:,j) = sqrt(nanmean((data.*areamat-refdata.*areamat).^2));
end

save(['RMSEarea_PFT_year_' num2str(year) '.mat'], 'RMSE','-v7.3');
