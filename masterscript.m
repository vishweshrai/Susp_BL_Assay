%   ***************** Quantification of SRB Release ***********************
%   ****** MASTER PROGRAM FOR RUNNING 'wellprops' AND GENERATING HIST *****
%   ****** OPTIMIZED FOR THE WELL SIZE OF 1 um  ***************************
%   _______________________________________________________________________
    clear
    clc
    close all
%   _______________________________________________________________________
%%  1. READING MASTER IMAGE
%   Directory of the images to be analyzed goes in 'dir'.
%   dir = "D:\Experimental_Data\Cytolysin_A\Suspended_Bilayer\SRB_Leakage_SULB_ClyA_Intermediates\Ternary_Lipids\SRB_PC_SM_Chol_1_um_well_ClyA_min_15";
%   Or use the following prompt to get the directory string:
%   dir = input("Input the directory where you've stored all the time-point images in one experiment:")
%   File name of the initial image goes in " ".
%   Name all the images in the sequence as per the time point, as the code does not have any segment to extract time-related meta data. :')
%   e.g., You can name the sequence as 0 min, 1 min, 5 min, etc.
%   Specify the first image in the sequence below:
%   first = input("Input filename of the first image (be careful of the case):")
%   first = "\" + first;
    I = imread(dir + first);
%   _______________________________________________________________________
%%  2. Making directories for storing results and processed images
    resdir = dir+"\Results";
    regimages = resdir+"\Registered_Images";
    regoutline = resdir+"\Registered_and_Outlined";
    if isfolder(resdir) == 1
       rmdir(resdir,"s");
    end
    mkdir(resdir);
    mkdir(regimages);
    mkdir(regoutline);
%   _______________________________________________________________________
%%  3. GENERATING STATISTICS FOR MASTER IMAGE
    prop = wellprop(I,resdir);
%   _______________________________________________________________________
%%  4. CREATING IMAGE DATASTORE FOR ALL IMAGES IN ONE RUN
    ds = imageDatastore(dir);
    fullFileNames = vertcat(ds.Files);
    [folder,FileName,ext] = fileparts(fullFileNames);
    mov = cell(numel(ds.Files),1);
    status = zeros(numel(ds.Files),1);
    lostarea = zeros(numel(ds.Files),1);
    regmask = cell(numel(ds.Files),1);
    medint = zeros(numel(ds.Files),1);
    time = zeros(numel(ds.Files),1);
    for i = 1:numel(ds.Files)        
        str = string(FileName(i));
        c1 = regexp(str,'[0-9]',"match");
        c2 = regexp(str,'sec',"match");
        if numel(c1) ~= 0
            s = join(c1,'');
            n = str2double(s);
        else
            n = 0;
        end
        if c2 == "sec"
           time(i) = n/60;
        else
           time(i) = n;
        end
    end
%   _______________________________________________________________________
%%  5. REGISTERING SUBSEQUENT IMAGES WITH MASTER IMAGE, AND MODIFYING
%   MASTER MASK TO ACCOUNT FOR MAX AREA LOST TO REGISTRATION.
    fix = I;
    fixbg = imgaussfilt(fix,50);
    fix = imsubtract(fix, fixbg);
    fix = imadjust(fix);
    panic = 0;
    route = 1;
    success = 0;
    Rfix = imref2d(size(fix));
    while (panic <= 1)&&(success == 0)
        for i = 1:numel(ds.Files)
            mov{i} = read(ds);
            movbg = imgaussfilt(mov{i},50);
            mov{i} = imsubtract(mov{i}, movbg);
            disp("Registering "+FileName(i)+" ...")
            movreg = wellreg(mov{i},fix,panic);
            mov{i} = movreg.moved;
            lost = lostfield(mov{i});
            lostarea(i,1) = lost.percent;
            regmask{i} = lost.mask;
        end
        if (max(lostarea) >= 0.30)&&(route == 1)
            panic = 1;
            route = route + 1;
            reset(ds);
            continue
        elseif (max(lostarea) >= 0.30)&&(route == 2)
            panic = 2;
        elseif (max(lostarea) <= 0.30)
            success = 1;
        end
    end
    if (max(lostarea) >= 0.30)&&(panic == 2)
        disp("Losing too much field of view, consider changing the" + ...
            "registration algorithm or revising match parameters!")
        return
    end
    for i = 1:numel(ds.Files)
        imwrite(mov{i},regimages+"\reg_"+FileName(i)+".tif",'tif');
        outline{i} = mov{i};
        imshow(outline{i},[]);
        hold on
        viscircles(prop.centres, prop.rad,"LineWidth",0.25,"Color","green")
%       outline{i} = insertShape(outline{i},"Circle",[prop.centres prop.rad], ...
%           "LineWidth",1,"Color","green",'Opacity',0.7);
%       imshow(outline{i},[]);
        hold off
        g(6) = gcf;
        exportgraphics(g(6),regoutline+"\regout_"+FileName(i)+".tif","Resolution",600);
        worstmask = ~(zeros(size(fix)));
    end
    for i = 1:numel(ds.Files)
        worstmask = and(worstmask,regmask{i});
    end
    figure('units','normalized','outerposition',[0 0 1 1])
    imshow(worstmask,[])
    mask = prop.mask;
    worstmask(~mask) = 0;
    figure('units','normalized','outerposition',[0 0 1 1])
    imshow(worstmask,[])
    wells = regionprops("table",worstmask,worstmask,"Circularity");
%   _______________________________________________________________________
%%  6. APPLYING MODIFIED MASK TO ALL IMAGES IN THE SERIES, AND EXTRACTING
%   INTENSITIES OF ALL REGIONS
    roll = 50;
    intensity = zeros(numel(ds.Files),numel(wells.Row));
    for k = 1:numel(ds.Files)
    %   Subtracting background intensities to estimate absolute fluorescent
    %   intensity 
    %         bg = imgaussfilt(mov{k},roll);
    %         mov{k} = imsubtract(mov{k},bg);
        mov{k}(~worstmask) = 1;
        stats = regionprops("table",worstmask,mov{k},"MeanIntensity");
        intensity(k,1:numel(stats.MeanIntensity)) = transpose( ...
            stats.MeanIntensity);
        medint(k) = median(intensity(k,:));
    end
    [time, order] = sort(time);
    intensity = intensity(order,:);
    medint = medint(order);
    %   Normalizing intensities wrt the median intensity at t = 0
    intensity(:,:) = intensity(:,:)/medint(1);
%   _______________________________________________________________________
%%  7. GENERATING INTENSITY VS TIME CURVES
    figure('units','normalized','outerposition',[0 0 1 1])
    for i=1:numel(intensity(1,:))
        plot(time,intensity(:,i));
        hold on
    end
    xlabel('Time [min]')
    ylabel('Mean Fluorescence Pixel Intensity of Well')
    title('Tracking Temporal Variations in Intensities of Single-Wells')
    g(1) = gcf;
    exportgraphics(g(1),resdir+"\single_well.jpg","Resolution",600);
    hold off
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(time,medint/medint(1),"LineWidth",2)
    xlabel('Time [min]')
    ylabel('Median Fluorescence Pixel Intensity Across All Wells')
    title('Temporal Variation of Median Intensity Across All Wells')
    g(2) = gcf;
    exportgraphics(g(2),resdir+"\average_well.jpg","Resolution",600);
    filename = resdir+"\Intensity_data.xlsx";
    writecell({'Time'},filename,'Sheet',1,'Range','B2');
    writecell({'Individual Pore Intensity'},filename,'Sheet',1,'Range' ...
        ,'C2');
    writecell({'Circularity'},filename,'Sheet',1,'Range','B3');
    writematrix(transpose(wells.Circularity),filename,'Sheet',1,'Range' ...
        ,'C3');
    writematrix([time, intensity],filename,'Sheet',1,'Range','B4');
    lgd = cell(numel(ds.Files),1);
    FileName = FileName(order);
    figure('units','normalized','outerposition',[0 0 1 1])
    for i = 1:numel(ds.Files)
        [heights, edges] = histcounts(intensity(i,:),"Normalization" ...
            ,"pdf");
        t = linspace(edges(1)-0.5*(edges(2)-edges(1)),edges(end)+ ...
            0.5*(edges(2) ...
            -edges(1)),numel(edges));
        dt = diff(t);
        Fvals = cumsum([0,heights.*dt]);
        F = spline(t, [0, Fvals, 0]);
        DF = fnder(F);
        fnplt(DF,'-',1)
        lgd{i} = string(FileName(i));
        hold on
    end
    hold off
    xlabel('Intensity');
    ylabel('Probability Density');
    title("Tracking Probability Densities of Well Intensities");
    legend(lgd)
    g(3) = gcf;
    exportgraphics(g(3),resdir+"\pdf.jpg","Resolution",600);
%   -----------------------------------------------------------------------
    nt = numel(intensity(:,1));
    ni = numel(intensity(1,:));
    edges = linspace(0,2,201);
    for j = 1:nt
        h(:,j) = histcounts(intensity(j,1:ni),edges,"Normalization","probability");
        h(:,j) = smooth(smooth(h(:,j),'moving'),'moving');
    end
    for k = 1:numel(edges)-1
        x(k) = edges(k)+0.5*(edges(k+1)-edges(k));
    end
    labels = num2cell(time(:));
    figure('units','normalized','outerposition',[0 0 1 1])
    y = 1:length(time(:));
    [X,Y] = meshgrid(x,y);
    w = waterfall(X,Y,transpose(h))
    set(w, 'FaceColor', 'flat');
    set(w, 'EdgeColor', 'k');
    w.EdgeAlpha = 0.8;        
    set(w, 'FaceVertexCData', rand(nt,3))
    w.FaceAlpha = 0.5;
    set(gca,'XGrid','on','YGrid','on','ZGrid','on')
    xlabel('Intensities')
    ylabel('Time [min]')
    yticks([1:(length(labels))])
    yticklabels(labels)
    zlabel('Probability Density')
    title(['Probability Distribution of Individual Well Intensities' ...
        ' with Time'])
    %pbaspect([2 4 1])
    g(4) = gcf;
    exportgraphics(g(4),resdir+"\3d_pdf.jpg","Resolution",600);
%   -----------------------------------------------------------------------
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc(1:numel(time),edges,h);
    Color = load("CustomColorMap.mat");
    ColorMap = Color.CustomColormap;
    colormap(ColorMap)
    colorbar
    xticks([1:(length(labels))])
    xticklabels(labels)
    labint = linspace(0,2,201);
    labint = num2cell(labint(:));
    xlabel('Time [min]');
    ylabel('Normalized Intensities');
%   title(['Probability Distribution of Individual Well Intensities' ...
%         ' with Time'])
    fontsize(gca,42,"pixels")
    set(gca,"FontWeight","bold")
    %daspect([3 50.0 1])
    set(gca,'Ydir','normal')
    g(5) = gcf;
    exportgraphics(g(5),resdir+"\2d_pdf.tif","Resolution",600);
%   _______________________________________________________________________
%
