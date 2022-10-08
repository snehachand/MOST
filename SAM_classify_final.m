clc;clearvars;close all;format compact; 
S1 = [];
S2 = [];
S3 = [];
S1_b = [];
S2_b = [];
S3_b = [];

%extracting black frame
fid_black = fopen('D:\Downloads\oral_database\reference\1_unprocessed\context\non_uniformity\dark_reference.raw','r'); 
data_black = fread(fid_black,'float');
fclose(fid_black);
img_black = reshape(data_black,2048,1088)';
black = zeros(272,512,16);
for i = 1:16
    black(:,:,i) = img_black(rem(i-1,4)+1:4:end,ceil(i/4):4:end);
end
black_vec = reshape(black, [size(black,1)*size(black,2),size(black,3)]);

image_list = dir('*.mat');
for i=1:1:length(image_list)
    load(image_list(i).name)
    multicube = Ms_cube;
    for l = 1:1:length(Mask_s)
        if isempty(Mask_s(l))==0
            Hist_g = reshape(multicube, [size(multicube,1)*size(multicube,2),size(multicube,3)]);
            black_g = Hist_g./black_vec;
            mask_ls = reshape (Mask_s{l}, [size(Mask_s{l},1)*size(Mask_s{l},2),1]);
            Hist_g = Hist_g(find(mask_ls==1),:);
            black_g = black_g(find(mask_ls==1),:);

            switch l
                case 1
                    S1 = [S1;Hist_g];
                    S1_b = [S1_b; black_g];
                case 2
                    S2 = [S2;Hist_g];
                    S2_b = [S2_b; black_g];
                case 3 
                    S3 = [S3;Hist_g];
                    S3_b = [S3_b;black_g];
                    
            end
        end
    end
end
%reflectance
ref_S1=S1./255;
ref_S2=S2./255;
ref_S3=S3./255;
%absorbance
absn_S1 = -1*log(ref_S1./S1_b);
absn_S2 = -1*log(ref_S2./S2_b);
absn_S3 = -1*log(ref_S3./S3_b);

S{1}=absn_S1;
S{2}=absn_S2;
S{3}=absn_S3;

WL_list=linspace(460, 600, 16);
 for j= 1:1:length(Mask_s)
    H_graph{j} = S{j};
    m_g{j}=mean(abs(H_graph{j}));
    std_g{j}=std(abs(H_graph{j}),0,1);
 end
%% SMA analysis
close all

%Import sample signal for which you want to plot SAM classification
%fid = fopen('E:\MATLAB\Oral_Database\raw_image\P1\unprocessed\img3.raw','r'); 
fid = fopen("D:\Downloads\oral_database\p5_unprocessed\p2\image_0000000000.raw",'r');
%fid = fopen("D:\Downloads\oral_database\p2_1_unprocessed\sathish\image_0000000000.raw",'r');
data = fread(fid,'float');
fclose(fid);
img = reshape(data,2048,1088)';

%extracting balck frame
fid_black = fopen('D:\Downloads\oral_database\reference\1_unprocessed\context\non_uniformity\dark_reference.raw','r'); 
data_black = fread(fid_black,'float');
fclose(fid_black);
img_black = reshape(data_black,2048,1088)';
black = zeros(272,512,16);
fst = zeros(272,512,16);
for i = 1:16
    fst(:,:,i) = img(rem(i-1,4)+1:4:end,ceil(i/4):4:end);
    black(:,:,i) = img_black(rem(i-1,4)+1:4:end,ceil(i/4):4:end);
end

%reflectance
ref_fst = fst./255;
%absorbance
absn_fst = -1*log(ref_fst./black);

%% Choose reference signal
%S ={[leukoplakia],[SCC],[normal]};

ref_signal_t=mean(abs(S2));
ref_signal_h=mean(abs(S3));
sample_signal = reshape(fst,[size(fst,1)*size(fst,2),size(fst,3)]);

%SAM angle : Add a relevant SAM angle, smaller the SAM angle, closer the
%sample and ref signal
SAM_Angle = 4;
%selected sample signal; t_sig = test signal or signal under test
t_sig=abs(sample_signal); 

 o1=0;
 m1=0;
 m2=0;
 for kk=1:1:size(t_sig,2)
        o1=o1+repmat(ref_signal_t(1,kk),[size(t_sig,1) 1]).*t_sig(:,kk);
        m1=m1+repmat(ref_signal_t(1,kk),[size(t_sig,1) 1]).^2;
        m2=m2+t_sig(:,kk).^2;
 end
 p1=0;
 q1=0;
 q2=0;
 for kk=1:1:size(t_sig,2)
        p1=p1+repmat(ref_signal_h(1,kk),[size(t_sig,1) 1]).*t_sig(:,kk);
        q1=q1+repmat(ref_signal_h(1,kk),[size(t_sig,1) 1]).^2;
        q2=q2+t_sig(:,kk).^2;
 end

%sam equation
sam_t=acos(o1./(m1.^0.5.*m2.^0.5)); 
sam_h=acos(p1./(q1.^0.5.*q2.^0.5)); 
%radian to degree
sam_t = (sam_t.*180)./3.14;
sam_h = (sam_h.*180)./3.14;

%creating rgb image of the data
rgb_img(:,:,1)=fst(:,:,1);
rgb_img(:,:,2)=fst(:,:,2);
rgb_img(:,:,3)=fst(:,:,4);
imshow(uint8(rgb_img));  

sam_imt_idx = find(sam_t<=SAM_Angle); %red
sam_imh_idx = find(sam_h<=SAM_Angle); %green
sam_img = uint8(rgb_img);
sam_img(sam_imt_idx) = 255; sam_img(sam_imt_idx+272*512) = 0; sam_img(sam_imt_idx+272*512*2) = 0;
sam_img(sam_imh_idx) = 0; sam_img(sam_imh_idx+272*512) = 255; sam_img(sam_imh_idx+272*512*2) = 0;
figure(2)
imshow(sam_img);

%% graph
absn_fst = reshape(absn_fst, 272*512,16);
t_a = absn_fst(find(sam_t<=SAM_Angle),:);
t_mean = mean(t_a);
t_std = std(t_a);

h_a = absn_fst(find(sam_h<=SAM_Angle),:);
h_mean = mean(h_a);
h_std = std(h_a);

    figure(3)
    x=linspace(460,600,16);
    plot(x,t_mean,'r',x, h_mean,'g', 'LineWidth', 4);
    hold on
    axis([460 600 3.4 5]);
    a=area(x,[t_mean'-t_std',2*t_std'])
    alpha(0.1)
    set(a, 'edgecolor','none')
    hold on 
    b=area(x,[h_mean'-h_std',2*h_std'])
    alpha(0.1)
    set(b, 'edgecolor','none')
    newcolor1=[1 1 1;0 1 0;1 1 1;1 0 0];
    colororder(newcolor1)
    legend('mean spectrum for SCC','mean spectrum for healthy oral tissue','fontweight','bold','fontsize',20,'FontName','Arial')
    box off
    set(gca,'FontSize',20, 'FontName','Arial', 'FontWeight','bold','LineWidth',3)
    ylabel('Absorbance','fontweight','bold','fontsize',20,'FontName','Arial')
    xlabel('Wavelength(nm)', 'fontweight','bold','fontsize',20,'FontName','Arial')
    hold on
