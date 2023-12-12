function imagePatch_test
im=imread('peppers.png');
figure(1)
imshow(im)
p=round(ginput(1));
disp(p)
figure(2)
imshow(imagePatch(im,p,64))
