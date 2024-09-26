%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: counties_cyclones.m
% By: Stephie Fried
% Date: Summer 2021
% Purpose: Computes the percent of counties that experience a cyclone that
% have a PDD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

bpath = 'C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\process_output';
datapath = 'C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\data\raw_data';
statapath = 'C:\Users\sdfried\Dropbox (ASU)\Research\Adaptation3\submission\data\output\excel\';

options = optimset('Algorithm','levenberg-marquardt','TolFun',10e-30,'TolX',10e-30, ...
    'MaxFunEvals', 10000, 'MaxIter', 10000, 'Display','none');

cd(bpath);

%% Import data

%Load list of counties that could experience cyclones
file2 = [statapath, 'possible_fips.xlsx'];
fips_vec = xlsread(file2);

%Load list of counties-cyclone pairs that declare pdds
file1= [statapath, 'pdds.xlsx'];
M1  = xlsread(file1);
fips_pdd = M1(:,1); 
stormid_pdd =M1(:, 2); 


%Load storm track data
file3 = [statapath, 'cyclone_wind_speeds_all.xlsx'];
M3 = xlsread(file3);
year_cyclone = M3(:,1);
month_cyclone = M3(:, 2);
day_cyclone = M3(:, 3);
time_cyclone = M3(:, 4);
msw_cyclone = M3(:, 5);
lat_cyclone = M3(:, 6);
lon_cyclone = -M3(:, 7);
storm_id_vec = M3(:, 20);

%Wind bands are measured in nautical miles
NE_34_raw = M3(:, 8);
SE_34_raw = M3(:, 9);
SW_34_raw = M3(:, 10);
NW_34_raw = M3(:, 11);
NE_50_raw = M3(:, 12);
SE_50_raw = M3(:, 13);
SW_50_raw = M3(:, 14);
NW_50_raw = M3(:, 15);
NE_64_raw = M3(:, 16);
SE_64_raw = M3(:, 17);
SW_64_raw = M3(:, 18);
NW_64_raw = M3(:, 19);
max_wind_radii = M3(:,21);

wb_vec = {"max", "64", "50", "34"};
direct_vec ={"NE", "NW", "SW", "SE"};

%Interpolate -999 entries
for wbi = 2:1:4
    wb = wb_vec{wbi};
    for directi = 1:4
        direct = direct_vec{directi};
        eval(strcat('direct_wb = ',direct, '_', wb, '_raw;'));
        for i = 2:1:length(NE_34_raw)-1
            if direct_wb(i) <0
                if direct_wb(i-1) >=0 && direct_wb(i+1) >=0
                    direct_wb(i) = mean([direct_wb(i-1), direct_wb(i+1)]);
                elseif direct_wb(i-1) >=0 && direct_wb(i+1) <=0
                    direct_wb(i) = direct_wb(i-1);
                elseif direct_wb(i-1) <=0 && direct_wb(i+1)>=0
                    direct_wb(i) = direct_wb(i+1);
                end
            end
            
            eval(strcat(direct, '_', wb, ' = direct_wb', ';'));
        end
    end
end

%Interpolate zero entries when other directions are non-zero
for wbi = 2:1:4
    wb = wb_vec{wbi};
    for i =1:length(NE_34)
        %create a matrix of wind radii at point i
        eval(strcat('temp = [SE_', wb, '(i), SW_', wb, '(i), NW_', wb, '(i), NE_', wb', '(i)];'))
        
        %If at least one radii is greater than zero, then check that all
        %radii are greater than zero
        if max(temp)>0
            for directi = 1:4
                direct = direct_vec{directi};
                eval(strcat('direct_wb = ',direct, '_', wb, ';'));
                %If one radii equals zero, replace it with the average of
                %the nonzero radii
                if direct_wb(i) ==0
                    direct_wb(i) = mean(temp(temp>0));
                    eval(strcat(direct, '_', wb, ' = direct_wb', ';'));
                end
            end
        end
    end
end

%replace missing entries for max_wind_radii with zeros
max_wind_radii(max_wind_radii <0) =0;
max_wind_radii(isnan(max_wind_radii)) =0;


%% Create wind bands for each storm track

%Bearings for wind bands in radians
bear_inc = pi/16;
bear_vec = [0:bear_inc:2*pi - bear_inc]';

%Max-wind (radius is the same in every direction)
lat_max = zeros(length(msw_cyclone), length(bear_vec) + 1);
lon_max = zeros(length(msw_cyclone), length(bear_vec) + 1);

for trackpti = 1:length(msw_cyclone)
    
    %Lat/lon at trackpt in degrees
    lat_deg = lat_cyclone(trackpti);
    lon_deg = lon_cyclone(trackpti);
    
    %Calculate coordinates of max wind band
    for i= 1:length(bear_vec)
        [lat_max(trackpti, i), lon_max(trackpti, i)] = endpoint_func(lat_deg, lon_deg, max_wind_radii(trackpti), bear_vec(i));
    end
    lat_max(trackpti, end) =lat_max(trackpti, 1); lon_max(trackpti, end) = lon_max(trackpti, 1);
    
end

%Other wind bands
lat_34 = zeros(length(msw_cyclone), length(bear_vec) + 1);
lon_34 = zeros(length(msw_cyclone), length(bear_vec) + 1);
lat_50 = zeros(length(msw_cyclone), length(bear_vec) + 1);
lon_50 = zeros(length(msw_cyclone), length(bear_vec) + 1);
lat_64 = zeros(length(msw_cyclone), length(bear_vec) + 1);
lon_64 = zeros(length(msw_cyclone), length(bear_vec) + 1);

tic
for trackpti = 1:length(msw_cyclone)
    
    lat_deg = lat_cyclone(trackpti);
    lon_deg = lon_cyclone(trackpti);
    
    for wbi = 2:1:4
        wb = wb_vec{wbi};
        
        eval(strcat('NE = NE_', wb, ';'));
        eval(strcat('NW = NW_', wb, ';'));
        eval(strcat('SE = SE_', wb, ';'));
        eval(strcat('SW = SW_', wb, ';'));
        
        radii_vec = zeros(size(bear_vec));
        
        %Impose wind radii at each given point.
        initial = find(bear_vec== pi/4);
        cons1 = length(bear_vec)/4;
        radii_vec(initial) = NE(trackpti);
        radii_vec(initial + cons1) =  SE(trackpti);
        radii_vec(initial + 2*cons1) = SW(trackpti);
        radii_vec(initial + 3*cons1) = NW(trackpti);
        
        %Start bisection
        radii_vec(1) = mean([radii_vec(initial), radii_vec(initial + 3*cons1)]);
        radii_vec(1+ cons1) = mean([radii_vec(initial), radii_vec(initial + cons1)]);
        radii_vec(1+ 2*cons1) = mean([radii_vec(initial + cons1), radii_vec(initial + 2*cons1)]);
        radii_vec(1+ 3*cons1) = mean([radii_vec(initial + 2*cons1), radii_vec(initial + 3*cons1)]);
        
        for i = 3:4:length(radii_vec)-2
            radii_vec(i) = mean([radii_vec(i + 2), radii_vec(i - 2)]);
        end
        radii_vec(end-1) = mean([radii_vec(end -3), radii_vec(1)]);
        
        for i = 2:2:length(radii_vec)-1
            radii_vec(i) = mean([radii_vec(i + 1), radii_vec(i - 1)]);
        end
        radii_vec(end) = mean([radii_vec(end -1), radii_vec(1)]);
        
        %Calculate coordinates of wimind bands.
        lat_wb = zeros(33, 1); lon_wb = zeros(33,1);
        for i= 1:length(radii_vec)
            [lat_wb(i), lon_wb(i)] = endpoint_func(lat_deg, lon_deg, radii_vec(i), bear_vec(i));
        end
        lat_wb(end) =lat_wb(1); lon_wb(end) = lon_wb(1);
        
        eval(strcat('lat_', wb, '(trackpti, :)=lat_wb;'));
        eval(strcat('lon_', wb, '(trackpti, :)=lon_wb;'));
        
    end
end
toc



%%
msw_mat = zeros(length(fips_vec), 70);
for fipsi = 1:length(fips_vec)
    fipsi
    fips = fips_vec(fipsi);
    cd(datapath);
    if fips < 10001
        county = shaperead('c_11au16.shp', 'UseGeoCoords', true, 'Selector', {@(FIPS)strcmpi(FIPS, strcat('0', num2str(fips))),'FIPS'});
    else
        county = shaperead('c_11au16.shp', 'UseGeoCoords', true, 'Selector', {@(FIPS)strcmpi(FIPS, num2str(fips)),'FIPS'});
    end
    
    if isempty(county)==0 %found fips code with county
        lat_county = county.Lat;
        lon_county = county.Lon;
        lat_county_cent = county.LAT;
        lon_county_cent = county.LON;
    end
    cd(bpath);
    
    for storm_id = 1:70
        trackpt1  = find(storm_id_vec == storm_id);
        %Find distance between county centroid and cyclone center
        dist = zeros(size(trackpt1));
        for j = 1:1:length(trackpt1)
            trackpti = trackpt1(j);
            lat_deg = lat_cyclone(trackpti);
            lon_deg = lon_cyclone(trackpti);
            dist(j) = distance(lat_deg, lon_deg, lat_county_cent, lon_county_cent); %measures distance in arclength
        end
        %Convert distance to nautical miles
        miles_dist=deg2nm((dist))*1.15078;
        hello = find(miles_dist < 700);
        trackpt2 = trackpt1(hello);
        msw_vec = zeros(max(length(trackpt2), 1), 1);
        
        if isempty(trackpt2) ==0
            for j = 1:1:length(trackpt2)
                
                trackpti = trackpt2(j);
                lat_deg = lat_cyclone(trackpti);
                lon_deg = lon_cyclone(trackpti);
                
                %Check if the center of the storm is in county
                [in1,on1] = inpolygon(lon_deg,lat_deg, lon_county,lat_county);
                
                %If cyclone track is inside county,
                if in1==1 || on1 ==1
                    msw_vec(j) = msw_cyclone(trackpti);
                else
                    
                    %If cyclone track is not inside county, check whether cyclone track is
                    %in each wind band
                    for wbi = 1:4
                        wb = wb_vec{wbi};
                        
                        eval(strcat('lat_wb = lat_', wb, '(trackpti, :);'));
                        eval(strcat('lon_wb = lon_', wb, '(trackpti, :);'));
                        
                        %Does county intersect wind band?
                        [xint,yint] = polyxpoly(lon_county,lat_county,lon_wb,lat_wb);
                        
                        %Is county fully inside wind band?
                        [in,on] = inpolygon(lon_county_cent,lat_county_cent, lon_wb,lat_wb);
                        
                        %If county intersects wind band or county is fully
                        %inside wind band, assign it the wind speed that
                        %corresponds to that wind band
                        if isempty(xint) ~= 1 || isempty(yint) ~=1 || in ==1 || on ==1
                            if wbi > 1
                                msw_vec(j) = max(msw_vec(j), eval(wb));
                            else
                                msw_vec(j) = max(msw_vec(j), msw_cyclone(trackpti));
                            end
                            break; %don't need to check the smaller wind bands
                        end
                        
                    end%end of wb loop
                end %end of cyclone in county if
                if msw_vec(j) >=64
                    break;
                end
                
            end %end of track point vec loop
        end
        
        msw_mat(fipsi, storm_id) = max(msw_vec);
    end
end
toc

%% Fraction of counties that experience a category 1 or higher cyclone that have a PDD 
strikes =0; 
pdds=0; 

for fipsi = 1:length(fips_vec)
    fips = fips_vec(fipsi);
    for storm_id = 1:max(storm_id_vec)
        if msw_mat(fipsi, storm_id) >= 64
            strikes = strikes +1; 
            temp = find(stormid_pdd ==storm_id); 
            if isempty(temp)==0
                temp2 = find(fips_pdd(temp)==fips); 
                if isempty(temp2)==0
                    pdds = pdds+1; 
                end
            end
        end
    end
end

pdds/strikes

save results_0121



            