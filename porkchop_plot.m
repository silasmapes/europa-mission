%% ======= Pork-Chop Plot ======= %%
% -------------------------------------------------------------------------
% A Pork-Chop Plot shows the departure and arrival dates of a heliocentric
% orbit transfer.
% -------------------------------------------------------------------------

clear; clc; close all
load rEar.mat; load vEar.mat; load rJup.mat; load vJup.mat;

depInt = 180;
arrInt = 1400;

jdDep = juliandate(2034,4,1) + (0:depInt-1);
jdArr = juliandate(2034,4,1) + (0:arrInt-1);

depDates = datetime(jdDep,'ConvertFrom','juliandate');
arrDates = datetime(jdArr,'ConvertFrom','juliandate');

dt = zeros(depInt,arrInt);
C3 = NaN(depInt,arrInt);
vinf = NaN(depInt,arrInt);

for l = 1:depInt
    for m = 1:arrInt
        t = jdArr(m)-jdDep(l);
        if t > 0
            dt(l,m) = (jdArr(m) - jdDep(l))*86400;
        end
    end
end



for n = 1:depInt
    for o = 1:arrInt
        if dt(n,o) > 0
            
            [v1s, v2s] = lambert_goodings(rEar(n,:),rJup(o,:),dt(n,o),132712440018,0,true,3);
            [v1l, v2l] = lambert_goodings(rEar(n,:),rJup(o,:),dt(n,o),132712440018,0,false,3);
            
            C3s = norm(v1s' - vEar(n,:))^2;
            C3l = norm(v1l' - vEar(n,:))^2;

            if C3s < C3l
                C3(n,o) = C3s;
                vinf(n,o) = norm(v2s);
            else
                C3(n,o) = C3l;
                vinf(n,o) = norm(v2l);
            end

        end
    end
end


%% ----- Plot Creator ----- %%

figure(); hold on
C3levels = [80 85 100 120 140 180 200 250 400 600];
[c1,h1] = contour(jdDep,jdArr,C3',C3levels);
clabel(c1,h1)

colormap("turbo")
cb = colorbar;
clim([C3levels(1) C3levels(end)]);

vinflevels = [6 6.5 7 8 9 10 12 14 16 18 20 25 30 35 40];
[c2,h2] = contour(jdDep,jdArr,vinf',vinflevels,'LineStyle','--');
clabel(c2,h2)

tof = jdArr' - jdDep;
[c3,h3] = contour(jdDep,jdArr,tof,400:200:1400,'LineStyle',':');
clabel(c3,h3)

%%



function [rEar, vEar, rJup, vJup] = orbitData(dep,arr)
    for y = 1:depInt
        [rEar(y,:), vEar(y,:)] = planetEphemeris(dep(y),'Sun','Earth');
            
    end
    for z = 1:arrInt
        [rJup(z,:), vJup(z,:)] = planetEphemeris(arr(z),'Sun','Jupiter');
    end
end