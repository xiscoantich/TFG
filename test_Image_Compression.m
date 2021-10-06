function test_Image_Compression
%This function loads an image and compresses it and compares the results
%data previously saved.
%
load('test.mat');
test_img=img;
test_img.B=rgb2gray(imread('Images/Happy-Test-Screen.png'));
test_read_image(test_img,test);
test_FTT2(test_img,test);
test_cut(test_img,test);
test_uncompress(test_img,test);
end

function test_read_image(test_img, test)
error = norm(double(test.B)-double(test_img.B));
if error<=10e-6
  cprintf('green', 'test read image  pass \n');
else
  cprintf('err', 'test read image error \n');
end
end

function test_FTT2(test_img,test)
error = norm(double(test.Bt)-double(test_img.Bt));
if error<=10e-6
  cprintf('green', 'test FTT2  pass \n');
else
  cprintf('err', 'test FTT2 error \n');
end
end

function test_cut(test_img,test)
error = norm(double(test.Atlow)-double(test_img.Atlow));
if error<=10e-6
  cprintf('green', 'test cut  pass \n');

else
  cprintf('err', 'test cut error \n');
end
end

function test_uncompress(test_img,test)
error = norm(double(test.Alow)-double(test_img.Alow));
if error<=10e-6
  cprintf('green', 'test cut  pass \n');

else
  cprintf('err', 'test cut error \n');
end
end

