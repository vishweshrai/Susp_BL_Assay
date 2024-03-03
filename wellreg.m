function movreg = wellreg(mov,fix,panic)
RF = imref2d(size(fix));
tempmov = imadjust(mov);
if panic == 0
    fixPt = detectSIFTFeatures(fix,"EdgeThreshold",10,"Sigma",1);
    movPt = detectSIFTFeatures(tempmov,"EdgeThreshold",10,"Sigma",1);
    [f1,vpts1] = extractFeatures(fix,fixPt);
    [f2,vpts2] = extractFeatures(tempmov,movPt);
    indexPairs = matchFeatures(f1,f2,"Method",'Exhaustive', ...
        "MatchThreshold",10,"MaxRatio",0.5,Unique=true);
    matchedPoints1 = vpts1(indexPairs(:,1));
    matchedPoints2 = vpts2(indexPairs(:,2));
    [tForm,inlierindex,status] =  estimateGeometricTransform2D( ...
        matchedPoints2.Location, ...
        matchedPoints1.Location,"rigid","MaxDistance",1.5,"Confidence" ...
        ,99,"MaxNumTrials",1000000);
    moved = imwarp(mov,tForm,"OutputView",RF);
    movreg.inlierindex = inlierindex;
    movreg.status = status;
elseif panic == 1
    Rmoving = imref2d(size(mov));
    tform = imregcorr(tempmov,Rmoving,fix,RF,"translation");
    moved = imwarp(mov,tform,"OutputView",RF);
end
movreg.moved = moved;
end