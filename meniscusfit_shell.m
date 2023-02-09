function [outputcoef, rankout, errmax] = meniscusfit_shell(xinput, yinput, tol)

errsave = 1;
coefsave = [];
coefsave_copy2 = [];
rank = -1;
rank_copy2 = -1;
for n = 3:20
        coefinit = [xinput(length(xinput),1), zeros(1,n), 0.2509*2.2];
        coef_tmp = lsqnonlin(@(coefs) meniscusfit(xinput,yinput,coefs),coefinit);
        coef_tmp(1) = coef_tmp(1) + yinput(1) - meniscusfit(xinput(1),0,coef_tmp);
        errmax = (max(abs(meniscusfit(xinput,yinput,coef_tmp))));
        if errmax < errsave
            errsave = errmax;
            coefsave = coef_tmp;
            rank = n;
        end
end
errsave = 1;
for m = 3:20
    l_surf = sqrt(0.072/1e3/9.8)*1e3;
    coefinit = zeros(1,m+1);
    coefinit(m-2:m+1) = coefsave(rank-1:rank+2);
    coef_tmp = lsqnonlin(@(coefs) Copy_2_of_meniscusfit(xinput,yinput,coefs),coefinit);
    part_polynomial = polyval(coef_tmp,xinput);
    part_exp = exp(-xinput/l_surf).*sqrt(1./xinput);
    multpart = part_exp.*part_polynomial;
    constadd = yinput(1) - multpart(1);
%     coefinit = [zstar4(length(zstar4),1), zeros(1,n-1), 0.2509];
%     x = lsqnonlin(@(coefs) Copy_of_meniscusfit(rstar4,zstar4(:,1),coefs),coefinit);
    errmax = (max(abs(Copy_2_of_meniscusfit(xinput,yinput,coef_tmp))));
    %[0.00327970848360634,-0.0206677371623448,0.177969218711422,0.984686375766926]
        if errmax < errsave
            errsave = errmax;
            coefsave_copy2 = coef_tmp;
            rank_copy2 = m;
        end
end
    part_polynomial = polyval(coefsave_copy2,xinput);
    part_exp = exp(-xinput/l_surf).*sqrt(1./xinput);
    multpart = part_exp.*part_polynomial;
    constadd = yinput(1) - multpart(1);
rankout = rank_copy2;
outputcoef(1) = constadd;
outputcoef(2:rankout+2) = coefsave_copy2;
errmax = errsave;
end
        