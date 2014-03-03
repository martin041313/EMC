%load rgb image
src = 'lena.jpg';
rgbI = imread(src);

%rgbI = double(rgbI)/255;
%convert to lab
labTransformation = makecform('srgb2lab');
labI = applycform(rgbI,labTransformation);

%seperate l,a,b
l = labI(:,:,1);
a = labI(:,:,2);
b = labI(:,:,3);

rgbIt = double(rgbI)/255;
%convert to lab
%labTransformation = makecform('srgb2lab');
labIt = applycform(rgbIt,labTransformation);

%seperate l,a,b
lt = labIt(:,:,1);
at = labIt(:,:,2);
bt = labIt(:,:,3);

figure, imshow(l) , title('l');
figure, imshow(lt) , title('lt');


%figure, imshow(a) , title('a');
%figure, imshow(b) , title('b');


