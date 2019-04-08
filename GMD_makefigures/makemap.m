function m = makemap(lat,lon,data,res)
m=NaN(round(180/res),round(360/res));
for i=1:length(lat)
    x=floor((90-lat(i))/res)+1;
    y=floor((lon(i)+180)/res)+1;
    m(x,y)=data(i);
end
