function lost = lostfield(I)
    I = im2double(I);
    I = imcomplement(I);
    I = I == 1;
    I = imopen(I,strel("rectangle",[4 4]));
    I = imclearborder(I,1);
    I = imcomplement(I);
    I = imclearborder(I,1);
    I = imclose(I,strel("rectangle",[64 64]));
    I = imcomplement(I);
    I1 = imopen(I,strel("rectangle",[256 1]));
    I2 = imopen(I,strel("rectangle",[1 256]));
    I = and(~I1,~I2);
    stat = regionprops(I,"Area");
    s = size(I);
    a = s(1)*s(2);
    lost.percent = (a-(stat.Area))/a;
    lost.mask = I;
end
