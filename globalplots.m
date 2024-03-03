%   ***************** Quantification of SRB Release ***********************
%   ***************************** Global Plots ****************************
%   _______________________________________________________________________
clear
clc
close all
%   _______________________________________________________________________
%%  1. READING RESULT GENERATED XCEL DATA

%   Add directories of resultant xcel sheets of 4 oligomers below in the
%   order 1, 15, 40 and 60
masterdir = "D:\Current Data\Suspended_Bilayer\SRB_Leakage_SULB_ClyA_Monomer\1_um";
dir(1) = importdata(masterdir+"\Binary_ClyA_Monomer\Region_01\Results\Intensity_Data.xlsx");
dir(2) = importdata(masterdir+"\Ternary_ClyA_Monomer\R_02\Results\Intensity_Data.xlsx");
%dir(3) = importdata(masterdir+"SRB_PC_SM_Chol_1_um_well_ClyA_min_40" + ...
%    "\Results\Intensity_Data.xlsx");
%dir(4) = importdata(masterdir+"SRB_PC_SM_Chol_1_um_well_ClyA_min_60" + ...
%    "\Results\Intensity_Data.xlsx");
system = 'Binary and Ternary';
%   _______________________________________________________________________
%%  2. PLOTTING THE DATA
%colorstring = 'mrgb';
colorstring = 'mg';
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:2
    nt = numel(dir(i).data(2:end,2));
    ni = numel(dir(i).data(2,2:end));
    intensity(1:nt,1:ni,i) = dir(i).data(2:end,2:end);
    time(1:nt,i) = dir(i).data(2:end,1);
    for j = 1:nt
        med(j,i) = median(intensity(j,1:ni,i));
        error(j,i) = std(intensity(j,1:ni,i))/med(1,i);
        % error(j,i) = std(intensity(j,1:ni,i))/sqrt(numel
        % (intensity(1:nt,1:ni,i)));
    end
    med(1:nt,i) = med(1:nt,i)/med(1,i);
    %       f = fit(time(1:nt,i),med(1:nt,i),'smoothingspline');
    %       xtime = linspace(time(1,i),time(nt,i),100*nt);
    %       plot(xtime,feval(f,xtime),"LineWidth",1,"Color",colorstring(i))
    errorbar(time(1:nt,i),med(1:nt,i),error(1:nt,i),"LineWidth",2, ...
        "Color",colorstring(i));
    %       errorbar(med(1:nt,i),error(1:nt,i));
    hold on
    clearvars xtime
end
hold off
xlabel("Time [min]");
ylabel("Median Intensity Across All Wells");
%legend('ClyA+DDM I/MG', 'ClyA+DDM SO/P', ...
%    'ClyA+DDM LO','ClyA+DDM Pore');
legend('ClyA with Binary Membrane','ClyA with Ternary Membrane')
set(gca,"FontSize", 24, "FontWeight","bold");
xlim([0 26]);
ylim([-1.09 1.09]);
title("Tracking Median Temporal Variation in "+system+" " + ...
    "Membrane System",FontSize=24,FontWeight="bold");
subtitle("  ");
g(1) = gcf;
exportgraphics(g(1),masterdir+system+"_Median_Time.jpeg", ...
    "Resolution",600);
%   _______________________________________________________________________