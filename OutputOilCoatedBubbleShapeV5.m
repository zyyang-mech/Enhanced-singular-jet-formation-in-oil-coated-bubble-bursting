function OutputOilCoatedBubbleShapeV5(SetRadius, SetOilfrac, filename, parametersarray)
%% OutputOilCoatedBubbleShape version 3
% Coded by Zhengyu Yang of Prof. Feng's FIT Group, MechSE UIUC. Last
% modified on Oct 24, 2021.
%
% Modification Aug 20: increased the data density.
% Modification Oct 24: changed output fitting function
% Modification Dec 8: changed meniscus_fit functions and added output for
% intersection.
% Modification Dec 19: allows higher oil fraction
%
% INPUT
% SetRadius: The radius of the air bubble. It refers to the horizontal
% radius by default.
% SetOilfrac: The oil fraction of the oil-coated bubble. It refers to the
% ratio of oil volume and the sphere volume defined by SetRadius by
% default.
% filname: The prefix of the output filename.
% parametersarray: The physical parameters and different meanings of
% SetRadius and SetOilfrac. It is put in as a cell. The default setting is
% 10 cSt - water setup.
%

Radius = SetRadius;
if nargin<=2
    filename = 'Default';
end

GammaWater = 0.072;
GammaOA = 0.01914;
GammaOW = 0.0409;
flagmode = 2; % 1 for volume, 2 for hor radius
flagoilf = 2; % 1 for real ratio, 2 for ratio calc by radius
Densoil = 930;
if nargin == 4
        GammaWater = parametersarray(1);
        GammaOA = parametersarray(2);
        GammaOW = parametersarray(3);
        flagmode = parametersarray(4);
        flagoilf = parametersarray(5);
        Densoil = parametersarray(6);
end

GammaOAW = GammaOA+GammaOW;
GammaCap = GammaOAW+GammaWater;

Bo = 10^3*9.8*Radius^2/GammaWater;
BoEq1 = Densoil*9.8*Radius^2/GammaOA;
BoEq2 = (10^3-Densoil)*9.8*Radius^2/GammaOW;
BoEq3 = 10^3*9.8*Radius^2/GammaOAW;
BoEq4 = Bo;
l_surf = sqrt(0.072/1e3/9.8)*1e3;

difrad = 1; %init
difof  = 1;

bubbleodeoptions = odeset('MaxStep',0.001); % Make odesteps smaller. apply to eqs 1 and 3
farawayodeoptions = odeset('MaxStep',0.05); % Make odesteps smaller. apply to eq 4
setfittingerr = 5e-4;

%%

R0set = Radius;
zconnect = 0.4;
sgndifrad0 = 0;
sgndifof0 = 0;
radstep = 0.1;
ofstep = 0.05;
while (difrad>1e-6 || (difof>1e-6 && ofstep>eps))
    zini_star = 0.01;
    r_ini_star = sqrt(2-zini_star)*sqrt(zini_star);
    drdz_ini_star = (1-zini_star)/sqrt(2-zini_star)/sqrt(zini_star);
    rarrayini = [ r_ini_star; drdz_ini_star  ];
    tspan = [zini_star zconnect];
    [zstar1,rstar1] = ode45(@(z,r) Young_Laplace(z,r,BoEq1), tspan, rarrayini, bubbleodeoptions);

    sizesolve = size(zstar1);

    rconnect = rstar1(sizesolve(1),1);
    drconnect = rstar1(sizesolve(1),2);%- cot theta_c = drcri
    dzconnect = 1/drconnect; 
    
%     if drconnect<0
%         error('Sorry, we cannot solve high oil fraction for now. Please choose a lower oil fraction.')
%     end

    rini2 = 0.01;
    nbvpmesh = 200;
    rmesh = linspace(rconnect,rini2,nbvpmesh);
    pp1 = spline([0;rstar1(:,1)],[0;zstar1]);
%     pp2 = spline([0;rstar1(:,1)],[0;1./rstar1(:,2)]);
%     Boz2ini = -2;
%     z02ini = 0;
%     solinit2 = bvpinit(rmesh,@(x) [ppval(pp1,x);ppval(pp2,x)],[Boz2ini;z02ini]);
% 
%     sols2 = bvp4c(@(r,z,param) YL_rtoz(r,z,BoEq2,param(1)+param(2)*BoEq2),...
%              @(za,zb,param) boundz20(za,zb,param,zconnect,dzconnect,rini2),...
%              solinit2);
    hmism = 1;
    stopcrit = odeset('Events',@events);
    stoppingbubbleodeoptions = odeset(stopcrit,bubbleodeoptions);
    
    % sign = -1: estimated curvature too small, Bobott should be increased,
    % or thickness should be decreased. Vice versa Note: it looks like the
    % opposite, but whatever we are using is right
    sgnthickness0 = 0;
    sgnthicknessc = 0;
    
    thicstep = 1;
    hguess = 0.88*(SetOilfrac)^0.66;
    while(abs(hmism)>1e-6 && thicstep>1e-8)
        Bobottguess = -2-BoEq2*hguess;
        ztor2 = ode45(@(z,r) Young_LaplacewBoz0(z,r,BoEq2,Bobottguess),...
            [zconnect, -5], [rconnect, drconnect], stoppingbubbleodeoptions);
        switch ztor2.ie 
            case 1
                hmism = 0;
                sgnthicknessc = 1;
            case 2
                hend = ztor2.x(end);
                rend = ztor2.y(1,end);
                Rcap2 = 2/(BoEq2*hend - Bobottguess);
                %From mathematica the r now should be R/sqrt(101)
                hmism = Rcap2/sqrt(101) - rend;
%                 hmism = 0;
                sgnthicknessc = -1;
            case 3
                hmism = 1;
                sgnthicknessc = 1;
        end
        if (sgnthicknessc*sgnthickness0 == -1)
            thicstep = thicstep*0.5;
        end
        sgnthickness0 = sgnthicknessc;
        if (sgnthicknessc < 0)
            hguess = hguess - thicstep;
        else
            hguess = hguess + thicstep;
        end

    end
        
    zcri_star_left = zconnect;
    zcri_star_right = 1.95;
    range_zcri = zcri_star_right-zcri_star_left;
    %zcri_star = 1.52599;
    warnStruct = warning('off','MATLAB:ode45:IntegrationTolNotMet'); %Opened for debugging
    while(range_zcri>eps)

        zcri_star = (zcri_star_left+zcri_star_right)/2;
        zini_star3 = zconnect;
        r_ini_star3 = rconnect;
        drdz_ini_star = drconnect;
        rarrayini = [ r_ini_star3; drdz_ini_star  ];
        tspan = [zini_star3 zcri_star];
    % 
    %     %
    %     % Solve 
    %     % (r'^2+1)-rr''  
    %     % ------------- = 2 + Bo z
    %     % r(r'^2+1)^3/2
    %     %
    % 
    % sols2.param(1) -> -2/Rcap2
    % sols2.param(2) -> hend
        param1_sols2 = -2/Rcap2;
        param2_sols2 = hend-rend^2/(Rcap2+sqrt(Rcap2^2-rend^2));
        Boz0 = (10^3-Densoil)*9.8*Radius^2/GammaOAW*param2_sols2 -2*GammaOA/GammaOAW ...
                  +param1_sols2*GammaOW/GammaOAW;

        [zstar3,rstar3] = ode45(@(z,r) Young_LaplacewBoz0(z,r,BoEq3,Boz0), tspan, rarrayini,bubbleodeoptions);


        sizesolve3 = size(zstar3);

        rcri3 = rstar3(sizesolve3(1),1);
        drcri3 = rstar3(sizesolve3(1),2);%- cot theta_c = drcri
        dzcri3 = 1/drcri3;  

        sinThetacri = sqrt(1/(1+drcri3^2));
        Rcap = rcri3/sinThetacri;
    %     Bozinf = 2*GammaCap/GammaWater/Rcap -2*GammaOA/GammaWater + ...
    %         param1_sols2*GammaOW/GammaWater ...
    %         + BoEq4*(10^3-Densoil)/10^3*param2_sols2;
    %     zinfref = Bozinf/BoEq4;
        zinfref = 2*GammaCap/10^3/9.8/Radius^2/Rcap - 2*GammaOA/10^3/9.8/Radius^2 + ...
            param1_sols2*GammaOW/10^3/9.8/Radius^2 ...
            +(10^3-Densoil)/10^3*param2_sols2;
        Bozinf = BoEq4*zinfref;

        rinf=15;
        [rstar4, zstar4] = ode45(@(r,z) YL_rtoz(r,z,BoEq4,Bozinf),[rcri3;rinf],[zcri_star;dzcri3],farawayodeoptions);

        sizesolu2 = size(rstar4);
        zfaroff = zstar4(sizesolu2(1),1);
        if (zfaroff>zinfref)
            zcri_star_right = zcri_star;
            range_zcri=range_zcri/2;
        elseif (zfaroff<zinfref)
            zcri_star_left = zcri_star;
            range_zcri=range_zcri/2;
        else
            range_zcri=0;
        end
    %     
%         fprintf('%f\n',zcri_star); % Previously used to check convergence
    %         
    end

    bottomplotx = linspace(0,r_ini_star,10);
    bottomploty = 1-sqrt(1-bottomplotx.^2);
    capplotx = linspace(rcri3,0,30);
    capploty = zcri_star - sqrt(Rcap^2-rcri3^2) + sqrt(Rcap^2-capplotx.^2);
    bottomplotx2 = linspace(0,rini2,10);
    R = -2/param1_sols2;
    bottomploty2 = param2_sols2+R-sqrt(R^2-bottomplotx2.^2);
    
    upperval = ppval(pp1,rmesh);
    if drconnect>0
        lowerspline = spline(ztor2.y(1,:),ztor2.x);
        lowerval = ppval(lowerspline,rmesh);
        oilshape = upperval-lowerval;
        oilvolintd = pi*(rmesh(1:nbvpmesh-1).^2-rmesh(2:nbvpmesh).^2).*(oilshape(1:nbvpmesh-1)+oilshape(2:nbvpmesh))/2;
        oilvoltot = sum(oilvolintd);
    else
        outerspline = spline(ztor2.x,ztor2.y(1,:));
        innerspline = spline(zstar1,rstar1(:,1));
        oilvolouter = sum(-pi/3*diff(ztor2.x).*(ztor2.y(1,1:end-1).^2+ztor2.y(1,2:end).^2+ztor2.y(1,1:end-1).*ztor2.y(1,2:end)));
        oilvolinner = sum(pi/3*diff(zstar1).*(rstar1(1:end-1,1).^2+rstar1(2:end,1).^2+rstar1(1:end-1,1).*rstar1(2:end,1)));
        oilvoltot = oilvolouter - oilvolinner;
    end

    bubblevolbottom = 1/3*pi*zini_star^2*(3-zini_star);
    bubblevollowintd = 1/3*pi*(zstar1(2:sizesolve(1))-zstar1(1:sizesolve(1)-1)).* ...
        (rstar1(2:sizesolve(1)).^2+rstar1(1:sizesolve(1)-1).^2+rstar1(1:sizesolve(1)-1).*rstar1(2:sizesolve(1)))';
    bubblevollower = sum(bubblevollowintd);

    bubblevolupintd = 1/3*pi*(zstar3(2:sizesolve3(1))-zstar3(1:sizesolve3(1)-1)).* ...
        (rstar3(2:sizesolve3(1)).^2+rstar3(1:sizesolve3(1)-1).^2+rstar3(1:sizesolve3(1)-1).*rstar3(2:sizesolve3(1)))';
    bubblevolupper = sum(bubblevolupintd);
    capheight = Rcap-sqrt(Rcap^2-rcri3^2);  
    bubblevolcap = 1/3*pi*capheight^2*(3*Rcap-capheight);

    bubblevoltot = bubblevolbottom + bubblevollower + bubblevolupper + bubblevolcap;

    maxr = max(max(rstar3(:,1)),max(ztor2.y(1,:)));
    bubbler = (bubblevoltot*3/4/pi)^(1/3);
%     oilfrac = oilvoltot/(oilvoltot+bubblevoltot);
    switch flagmode
        case 1
            raddiml = bubbler;
            difrad = abs(Radius*bubbler-SetRadius);
            sgndifradc = sign(SetRadius-Radius*bubbler);
        case 2
            raddiml = maxr;
            difrad = abs(Radius*maxr-SetRadius);
            sgndifradc = sign(SetRadius-Radius*maxr);
    end
    switch flagoilf
        case 1
            oilfrac = oilvoltot/(oilvoltot+bubblevoltot);
            difof = abs(oilfrac-SetOilfrac);
        case 2
            oilfrac = oilvoltot/(4*pi/3*raddiml^3);
            difof = abs(oilfrac-SetOilfrac);
    end
    sgndifofc = sign(SetOilfrac-oilfrac);
    if (sgndifrad0*sgndifradc == -1)
        radstep = radstep*0.5;
    end
    if (sgndifof0*sgndifofc == -1)
        ofstep = ofstep*0.5;
    end
    sgndifrad0 = sgndifradc;
    sgndifof0  = sgndifofc ;
    if (sgndifradc>0)
        Radius = Radius+radstep;
    elseif (sgndifradc<0)
        Radius = Radius-radstep;
    end
    if (sgndifofc>0)
        zconnect = zconnect+ofstep;
    elseif (sgndifofc<0)
        zconnect = zconnect-ofstep;
    end
end
if abs(difof) > 1e-6
    warning(sprintf('difof is %f%% > 1e-6',abs(difof)));
end
warning(warnStruct);

filename_OW = strcat(filename,'_OW.csv');
filename_OA = strcat(filename,'_OA.csv');
filename_WA = strcat(filename,'_WA.csv');
filename_Inters = strcat(filename,'_intersection.csv');
f_OW = fopen(filename_OW,'w');
fprintf(f_OW,'%s\n','r/mm, z/mm');
fprintf(f_OW,'%s%.15f\n','0, ',param2_sols2*Radius*1e3);
f_OW = fclose(f_OW);
output_Eq2 = [    fliplr(ztor2.y(1,:))' fliplr(ztor2.x)']*Radius*1e3;
output_Eq3 = [rstar3(2:sizesolve3(1),1) zstar3(2:sizesolve3(1))]*Radius*1e3;
output_Eq3_real = [rstar3(1:sizesolve3(1),1) zstar3(1:sizesolve3(1))]*Radius*1e3; %The Eq3 that includes the section point
writematrix(output_Eq2, filename_OW,'WriteMode', 'append');
writematrix(output_Eq3, filename_OW,'WriteMode', 'append');
f_OA = fopen(filename_OA,'w');
fprintf(f_OA,'%s\n','r/mm, z/mm');
f_OA = fclose(f_OA);
output_Eq1 = [bottomplotx(1:9)' bottomploty(1:9)'
              rstar1(:,1)       zstar1]*Radius*1e3;
writematrix(output_Eq1, filename_OA,'WriteMode', 'append');
f_WA = fopen(filename_WA,'w');
fprintf(f_WA,'%s\n','r/mm, z/mm');
f_WA = fclose(f_WA);
output_Eq4 = [rstar4 zstar4(:,1)]*Radius*1e3;
writematrix(output_Eq4, filename_WA,'WriteMode', 'append');


center = [0, SetRadius]*1e3;
diff_Eq1 = output_Eq1-center;
theta_Eq1 = atan(diff_Eq1(:,1)./(-diff_Eq1(:,2)));
r_Eq1 = sqrt(diff_Eq1(:,1).^2+diff_Eq1(:,2).^2);
diff_Eq2 = output_Eq2-center;
theta_Eq2 = atan(diff_Eq2(:,1)./(-diff_Eq2(:,2)));
r_Eq2 = sqrt(diff_Eq2(:,1).^2+diff_Eq2(:,2).^2);
thetaend = theta_Eq1(end);
thetareg = linspace(0,thetaend,400);
rreg_Eq1 = spline(theta_Eq1,r_Eq1,thetareg);
rreg_Eq2 = spline(theta_Eq2,r_Eq2,thetareg);
thickreg = rreg_Eq2 - rreg_Eq1;

%     figure(2)
%     plot(rstar1(:,1),zstar1,'k');
%         hold on;
%     plot(ztor2.y(1,:),ztor2.x,'r');
%     plot(rstar3(:,1),zstar3,'k');
%     plot(rstar4,zstar4(:,1),'k');
%     plot(bottomplotx,bottomploty,'k');
%     plot(capplotx,capploty,'k');
%     plot(-rstar1(:,1),zstar1,'k');
%     plot(-ztor2.y(1,:),ztor2.x,'r');
%     plot(-rstar3(:,1),zstar3,'k');
%     plot(-rstar4,zstar4(:,1),'k');
%     plot(-bottomplotx,bottomploty,'k');
%     plot(-capplotx,capploty,'k');
%     axis equal;
    
    % Getting the fit polynomials
err = 1;
rankset = 7;
while (err>setfittingerr)
    rankset = rankset + 1;
    polyEq3 = polyalign(output_Eq3_real(:,2),output_Eq3_real(:,1),rankset,2);
    err = max(abs(output_Eq3_real(:,1) - polyval(polyEq3,output_Eq3_real(:,2))));
%     display(err);
end
filename_polyEq3 = strcat(filename,'_polyEq3.csv');
filename_polyEq3text = strcat(filename,'_polyEq3.txt');
writematrix(polyEq3, filename_polyEq3);
ft_polyEq3 = fopen(filename_polyEq3text, 'w');
if rankset > 99
    warning('The polynomial may have not been well printed');
end
for i = 1:rankset+1
    xpow = rankset+1-i;
    if xpow > 1
        fprintf(ft_polyEq3,'%+.15E%s%2i%s',polyEq3(i),'*pow(z,',xpow,'.)');
    elseif xpow == 1
        fprintf(ft_polyEq3,'%+.15E%s',polyEq3(i),'*z');
    elseif xpow == 0
        fprintf(ft_polyEq3,'%+.15E',polyEq3(i));
    end
end
fclose(ft_polyEq3);

% divide Eq1 to output
rInter2 = output_Eq3(end,1);
breakind = output_Eq1(:,1)<=rInter2;
theotherpart = ~(breakind&circshift(breakind,-1));
output_Eq1_part1 = output_Eq1(breakind,:);
output_Eq1_part2 = output_Eq1(theotherpart,:);

% Get polynomial for Eq1 part 1
err = 1;
rankset = 7;
while (err>setfittingerr)
    rankset = rankset + 1;
    polyEq1part1 = polyalign(output_Eq1_part1(:,1),output_Eq1_part1(:,2),rankset,2);
    err = max(abs(output_Eq1_part1(:,2) - polyval(polyEq1part1,output_Eq1_part1(:,1))));
end
filename_polyEq1part1 = strcat(filename,'_polyEq1part1.csv');
filename_polyEq1part1text = strcat(filename,'_polyEq1part1.txt');
writematrix(polyEq1part1, filename_polyEq1part1);
ft_polyEq1part1 = fopen(filename_polyEq1part1text, 'w');
if rankset > 99
    warning('The polynomial may have not been well printed');
end
for i = 1:rankset+1
    xpow = rankset+1-i;
    if xpow > 1
        fprintf(ft_polyEq1part1,'%+.15E%s%2i%s',polyEq1part1(i),'*pow(r,',xpow,'.)');
    elseif xpow == 1
        fprintf(ft_polyEq1part1,'%+.15E%s',polyEq1part1(i),'*r');
    elseif xpow == 0
        fprintf(ft_polyEq1part1,'%+.15E',polyEq1part1(i));
    end
end
fclose(ft_polyEq1part1);

% Get intersectionpoint for Eq.1
zInter0_Eq1 = polyval(polyEq1part1,rInter2);
output_Eq1_part2 = [[rInter2,zInter0_Eq1];output_Eq1_part2];

% Polynomial for Eq.1 part 2
err = 1;
rankset = 7;
while (err>setfittingerr)
    rankset = rankset + 1;
    polyEq1part2 = polyalign(output_Eq1_part2(:,2),output_Eq1_part2(:,1),rankset,2);
    err = max(abs(output_Eq1_part2(:,1) - polyval(polyEq1part2,output_Eq1_part2(:,2))));
end
filename_polyEq1part2 = strcat(filename,'_polyEq1part2.csv');
filename_polyEq1part2text = strcat(filename,'_polyEq1part2.txt');
writematrix(polyEq1part2, filename_polyEq1part2);
ft_polyEq1part2 = fopen(filename_polyEq1part2text, 'w');
if rankset > 99
    warning('The polynomial may have not been well printed');
end
for i = 1:rankset+1
    xpow = rankset+1-i;
    if xpow > 1
        fprintf(ft_polyEq1part2,'%+.15E%s%2i%s',polyEq1part2(i),'*pow(z,',xpow,'.)');
    elseif xpow == 1
        fprintf(ft_polyEq1part2,'%+.15E%s',polyEq1part2(i),'*z');
    elseif xpow == 0
        fprintf(ft_polyEq1part2,'%+.15E',polyEq1part2(i));
    end
end
fclose(ft_polyEq1part2);

% Repeat for Eq2
breakind2 = output_Eq2(:,1)<=rInter2;
theotherpart2 = ~(breakind2&circshift(breakind2,-1));
output_Eq2_part1 = output_Eq2(breakind2,:);
output_Eq2_part2 = output_Eq2(theotherpart2,:);

% Get polynomial for Eq2 part 1
err = 1;
rankset = 7;
while (err>setfittingerr)
    rankset = rankset + 1;
    polyEq2part1 = polyalign(output_Eq2_part1(:,1),output_Eq2_part1(:,2),rankset,2);
    err = max(abs(output_Eq2_part1(:,2) - polyval(polyEq2part1,output_Eq2_part1(:,1))));
end
filename_polyEq2part1 = strcat(filename,'_polyEq2part1.csv');
filename_polyEq2part1text = strcat(filename,'_polyEq2part1.txt');
writematrix(polyEq2part1, filename_polyEq2part1);
ft_polyEq2part1 = fopen(filename_polyEq2part1text, 'w');
if rankset > 99
    warning('The polynomial may have not been well printed');
end
for i = 1:rankset+1
    xpow = rankset+1-i;
    if xpow > 1
        fprintf(ft_polyEq2part1,'%+.15E%s%2i%s',polyEq2part1(i),'*pow(r,',xpow,'.)');
    elseif xpow == 1
        fprintf(ft_polyEq2part1,'%+.15E%s',polyEq2part1(i),'*r');
    elseif xpow == 0
        fprintf(ft_polyEq2part1,'%+.15E',polyEq2part1(i));
    end
end
fclose(ft_polyEq2part1);

% Get intersectionpoint for Eq.1
zInter0_Eq2 = polyval(polyEq2part1,rInter2);
output_Eq2_part2 = [[rInter2,zInter0_Eq2];output_Eq2_part2];

% Polynomial for Eq.2 part 2
err = 1;
rankset = 7;
while (err>setfittingerr)
    rankset = rankset + 1;
    polyEq2part2 = polyalign(output_Eq2_part2(:,2),output_Eq2_part2(:,1),rankset,2);
    err = max(abs(output_Eq2_part2(:,1) - polyval(polyEq2part2,output_Eq2_part2(:,2))));
end
filename_polyEq2part2 = strcat(filename,'_polyEq2part2.csv');
filename_polyEq2part2text = strcat(filename,'_polyEq2part2.txt');
writematrix(polyEq2part2, filename_polyEq2part2);
ft_polyEq2part2 = fopen(filename_polyEq2part2text, 'w');
if rankset > 99
    warning('The polynomial may have not been well printed');
end
for i = 1:rankset+1
    xpow = rankset+1-i;
    if xpow > 1
        fprintf(ft_polyEq2part2,'%+.15E%s%2i%s',polyEq2part2(i),'*pow(z,',xpow,'.)');
    elseif xpow == 1
        fprintf(ft_polyEq2part2,'%+.15E%s',polyEq2part2(i),'*z');
    elseif xpow == 0
        fprintf(ft_polyEq2part2,'%+.15E',polyEq2part2(i));
    end
end
fclose(ft_polyEq2part2);

[coefEq4, rankexpo, errexpo] = meniscusfit_shell(output_Eq4(:,1), output_Eq4(:,2));
tmpa=fopen('New.txt','a');
fprintf(tmpa,'%f\n',errexpo);
fclose(tmpa);
filename_coefEq4 = strcat(filename,'_coefEq4.csv');
filename_coefEq4text = strcat(filename,'_coefEq4.txt');
writematrix(coefEq4, filename_coefEq4);
ft_coefEq4 = fopen(filename_coefEq4text, 'w');
for i = 1:rankexpo + 2
    if i==1
        fprintf(ft_coefEq4,'%+.15E%s%+.15E%s',coefEq4(i),'+pow(2.718281828459046,-r/(',l_surf,'))*pow(1./r,0.5)*(');
    else
        xpow = rankexpo + 2 - i;
        if xpow > 1
            fprintf(ft_coefEq4,'%+.15E%s%2i%s',coefEq4(i),'*pow(r,',xpow,'.)');
        elseif xpow == 1
            fprintf(ft_coefEq4,'%+.15E%s',coefEq4(i),'*r');
        elseif xpow == 0
            fprintf(ft_coefEq4,'%+.15E%s',coefEq4(i),')');
        end
    end
end
fclose(ft_coefEq4);

f_ints = fopen(filename_Inters,'w');
fprintf(f_ints,'%s\n','Intersection equations,r/mm, z/mm');
fprintf(f_ints,'Eq.1_Inter,%.16f,%.16f\n',rInter2,zInter0_Eq1);
fprintf(f_ints,'Eq.2_Inter,%.16f,%.16f\n',rInter2,zInter0_Eq2);
fprintf(f_ints,'Eq.1-2-3,%.16f,%.16f\n',output_Eq1(end,1),output_Eq1(end,2));
if (output_Eq1(end,1)~=output_Eq2(end,1))
    error('Check ends meet')
end
fprintf(f_ints,'Eq.3-4,%.16f,%.16f\n',output_Eq3(end,1),output_Eq3(end,2));
if (output_Eq3(end,1)~=output_Eq4(1,1))
    error('Check ends meet (2)')
end
f_ints = fclose(f_ints);


% text(rinf-3,-0.75,sprintf('%s%d%s','Air volume ',bubblevoltot*Radius^3*1e9,' uL'));
% text(rinf-3,-1,sprintf('%s%d%s','Oil volume ',oilvoltot*Radius^3*1e9,' uL'));
% text(rinf-3,-0.25,sprintf('%s%d%s','Volumetric radius ',(bubblevoltot*3/4/pi)^(1/3)*Radius*1e3,'mm'));
% text(rinf-3,-0.5,sprintf('%s%d%s','Horizontal radius ',maxr*Radius*1e3,'mm'));
% text(rinf-3,0,sprintf('%s%.2f%s','Oil/all fraction ',oilfrac*100,'%'));
end

function res=Young_Laplace(z,r,Bofn)
res = [r(2);
       (-Bofn*z-2)*sqrt(r(2)^2+1)^3+((r(2)^2)+1)/r(1)];
end

function res=YL_rtoz(r,z,Bofn,Bozinf)
res = [z(2);
       sqrt(1+z(2)^2)^3*(Bofn*z(1)-Bozinf)-z(2)*(1+z(2)^2)/r];
end
function res=boundz20(za,zb,param,zacri1,dzacri1,rini)
R = -2/param(1);
res = [za(1)-zacri1;
       za(2)-dzacri1;
       zb(1)-param(2)-rini^2/(R+sqrt(R^2-rini^2));
       zb(2)-rini/sqrt(R^2-rini^2)];
end

function res=Young_LaplacewBoz0(z,r,Bofn,Boz0)
res = [r(2);
       (-Bofn*z+Boz0)*sqrt(r(2)^2+1)^3+((r(2)^2)+1)/r(1)];
end

function [value,isterminal,direction] = events(~,y)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = [y(1)-0.00001; y(2)-50; y(2)];     % detect r = 0 or large dr/dz
isterminal = [1; 1; 1];   % stop the integration
direction = [0;0;-1];   % negative direction
end