function prop = wellprop(Im,dir)
%   **** Function for Detecting Membrane Pores and Generating Binary Mask *
%   ***** TASK: Use Circular Hough Transform, make a mask containing ******
%   ***************** circles, apply mask on original image ***************
%   _______________________________________________________________________
%%   1. IMAGE PREPROCESSING
%   Resizing so that circular hough transform can be safely applied at ~5 
%   pixels. Adding more pixels to the image to effect detection of smaller 
%   circles. This is required ONLY for analyzing 1 micron wells.
%   Default intensity-interpolation route for added pixels is "bicubic".
%   Caveat: Resolution may go down due to interpolation.     
    Im = imresize(Im,2,"bicubic");
    original = Im;
    IBG = imgaussfilt(Im,100);
    Im = imsubtract(Im,IBG);
    IConAdj = imadjust(Im);
    IConAdj = adapthisteq(IConAdj);  %  Scrambling pixel-intensities to fit
                                                    %  uniform distribution
%   _______________________________________________________________________
%%   2. PRELIMINARY BINARIZATION
%   Automatic (custom-settings based) Binarization
    IBW = imbinarize(IConAdj,"adaptive","ForegroundPolarity", ...
        "bright","Sensitivity",0.15);
%   _______________________________________________________________________
%%   3. FAST MARCHING
    WGDiff = graydiffweight(IConAdj,IBW);
    IFastGDiff = imsegfmm(WGDiff,IBW,0.00001);
%   _______________________________________________________________________
%%   4. CLEANING BINARY MASK
%   Knocking-off regions with less-than-or-equal-to 5 connected pixels
    IFastGDiff = bwareaopen(IFastGDiff,5);
    figure('units','normalized','outerposition',[0 0 1 1])
    imshowpair(IFastGDiff,original);
    prelim = gcf;
    exportgraphics(prelim,dir+"\preliminary_binarization_master_" + ...
        "image.tif","Resolution",600);
%   _______________________________________________________________________
%%   5. FINDING CIRCLES IN THE BINARY IMAGE
    rmin = 9;
    rmax = 12;
    [c,r] = imfindcircles(IFastGDiff,[rmin rmax],"Method","TwoStage", ...
        "EdgeThreshold",0.05,"ObjectPolarity","bright","Sensitivity",0.9);
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(original,[])
    hold on
    viscircles(c,r,"LineWidth",0.25)
    hold off
%   _______________________________________________________________________
%%   6. GENERATING MASK USING DETECTED CIRCLES
    Mask = zeros(size(Im));
    Mask = insertShape(Mask,"FilledCircle",[c(:,:) r(:)]);
    Mask = im2gray(Mask);
    Mask = imbinarize(Mask);
    Mask = imopen(Mask,strel("diamond",3));
    cc = bwconncomp(Mask);
    original(~Mask) = 1;
    h2 = figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(original,[])
%   _______________________________________________________________________
%%   7. STORING EXTRACTED FEATURES AS FIELDS OF AN OUTPUT STRUCTURE
    prop.mask = imresize(Mask,0.5);
    prop.outline = h1;
    prop.masked = h2;
    prop.centres = c*0.5;
    prop.rad = r*0.5;
    prop.NumHoles = cc.NumObjects;
    prop.id = cc.PixelIdxList;
end