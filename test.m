I=imread('pic2.png');
I=double(I)/255;
%colorTransform = makecform('srgb2lab');
%lab = applycform(rgbImage, colorTransform);
w     = 5;       % bilateral filter half-width
sigma = [3 0.1]; % bilateral filter standard deviations

I1=bfilter2(I,w,sigma);

subplot(1,2,1);
imshow(I);
subplot(1,2,2);
imshow(I1)