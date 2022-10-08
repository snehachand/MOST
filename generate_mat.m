clc;clearvars;close all;format compact; %aw
%extracting raw image
fid = fopen("D:\oral_database\Cancer_Segmentation\raw_files\p5_unprocessed\p2\image_0000000000.raw",'r');

rawdata = fread(fid,'float');
fclose(fid);
img = reshape(rawdata,2048,1088)';
figure;imshow(img,[]);

Ms_cube = zeros(272,512,16);
for i = 1:16
    Ms_cube(:,:,i) = img(rem(i-1,4)+1:4:end,ceil(i/4):4:end);
end

%Mask_s ={[leukoplakia],[SCC],[normal]};
Mask_s ={[],[],[]};
maskedImage_s ={[],[],[]};
imageSegmenter( Ms_cube(:,:,1));
Mask_s{2}=BW;
maskedImage_s{2}=maskedImage;

figure;imshow(Mask_s{2},[])
figure;imshow(maskedImage_s{2},[])
imwrite(Ms_cube, 'test.tiff', 'WriteMode', 'append');
save('IMG.mat','rawdata','img','Ms_cube','Mask_s','maskedImage_s');
