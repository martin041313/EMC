I=imread('pic2.png');
%I=double(I)/255;
%colorTransform = makecform('srgb2lab');
%lab = applycform(rgbImage, colorTransform);


I1=emc(I,25,25);

subplot(1,2,1);
imshow(I);
subplot(1,2,2);
imshow(I1)