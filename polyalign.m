function res = polyalign(xinput,yinput,rank,num,boundaries)
%% Polyfit with boundary aligned to data boundary
size = length(xinput);
size_y = length(yinput);
if (size ~= size_y)
    error('vector size does not match');
end
res = polyfit(xinput, yinput, rank);
switch num
    case 2
        ydif = yinput([1,size]) - polyval(res,xinput([1,size]));
        pdif = polyfit(xinput([1,size]), ydif, 1);
        res(rank:rank+1) = res(rank:rank+1) + pdif;
    otherwise
end
end