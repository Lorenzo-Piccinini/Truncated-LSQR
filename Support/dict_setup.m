function[H1,H2,I1,I2,rhs1,rhs2]=dict_setup(n_sample,dim_patch,digit)

% dict_setup(n_sample,dim_patch)
% dict_setup(n_sample,n_pix,n_pix_rhs)
% H1 = zeros(n_sample*28, n_pix*10);
% H2 = zeros(n_sample*28, n_pix*10);
% I1 = zeros(n_sample*28, n_pix*10);
% I2 = zeros(n_sample*28, n_pix*10);
% 
% rhs1 = zeros(n_sample*28, n_pix_rhs);
% 
% addpath('./Data')
% load mnist_all
% 
% for k=1:n_sample
% 
%     im0 = double(train0(1+k,:)');
%     im0 = reshape(im0,28,28);
%     im1 = double(train1(1+k,:)');
%     im1 = reshape(im1,28,28);
%     im2 = double(train2(1+k,:)');
%     im2 = reshape(im2,28,28);
%     im3 = double(train3(1+k,:)');
%     im3 = reshape(im3,28,28);
%     im4 = double(train4(1+k,:)');
%     im4 = reshape(im4,28,28);
%     im5 = double(train5(1+k,:)');
%     im5 = reshape(im5,28,28);
%     im6 = double(train6(1+k,:)');
%     im6 = reshape(im6,28,28);
%     im7 = double(train7(1+k,:)');
%     im7 = reshape(im7,28,28);
%     im8 = double(train8(1+k,:)');
%     im8 = reshape(im8,28,28);
%     im9 = double(train9(1+k,:)');
%     im9 = reshape(im9,28,28);
% 
%     for j=1:28
% 
%         pix = randperm(28,n_pix);
% 
%         H1(j+(k-1)*28,:)=[im0(j,pix), im1(j,pix), im2(j,pix), im3(j,pix), im4(j,pix), im5(j,pix), im6(j,pix), im7(j,pix), im8(j,pix), im9(j,pix)];
%     end
% 
% 
% end
% 
% 
% for k=1:n_sample
% 
%     im0 = double(train0(400+k,:)');
%     im0 = reshape(im0,28,28);
%     im1 = double(train1(400+k,:)');
%     im1 = reshape(im1,28,28);
%     im2 = double(train2(400+k,:)');
%     im2 = reshape(im2,28,28);
%     im3 = double(train3(400+k,:)');
%     im3 = reshape(im3,28,28);
%     im4 = double(train4(400+k,:)');
%     im4 = reshape(im4,28,28);
%     im5 = double(train5(400+k,:)');
%     im5 = reshape(im5,28,28);
%     im6 = double(train6(400+k,:)');
%     im6 = reshape(im6,28,28);
%     im7 = double(train7(400+k,:)');
%     im7 = reshape(im7,28,28);
%     im8 = double(train8(400+k,:)');
%     im8 = reshape(im8,28,28);
%     im9 = double(train9(400+k,:)');
%     im9 = reshape(im9,28,28);
% 
%     for j=1:28
% 
%         pix = randperm(28,n_pix);
% 
%         H2(j+(k-1)*28,:)=[im0(j,pix), im1(j,pix), im2(j,pix), im3(j,pix), im4(j,pix), im5(j,pix), im6(j,pix), im7(j,pix), im8(j,pix), im9(j,pix)];
% 
%     end
% 
% end
% 
% for k=1:n_sample
% 
%     im0 = double(train0(800+k,:)');
%     im0 = reshape(im0,28,28);
%     im1 = double(train1(800+k,:)');
%     im1 = reshape(im1,28,28);
%     im2 = double(train2(800+k,:)');
%     im2 = reshape(im2,28,28);
%     im3 = double(train3(800+k,:)');
%     im3 = reshape(im3,28,28);
%     im4 = double(train4(800+k,:)');
%     im4 = reshape(im4,28,28);
%     im5 = double(train5(800+k,:)');
%     im5 = reshape(im5,28,28);
%     im6 = double(train6(800+k,:)');
%     im6 = reshape(im6,28,28);
%     im7 = double(train7(800+k,:)');
%     im7 = reshape(im7,28,28);
%     im8 = double(train8(800+k,:)');
%     im8 = reshape(im8,28,28);
%     im9 = double(train9(800+k,:)');
%     im9 = reshape(im9,28,28);
% 
%     for j=1:28
% 
%         pix = randperm(28,n_pix);
% 
%         I1(j+(k-1)*28,:)=[im0(j,pix), im1(j,pix), im2(j,pix), im3(j,pix), im4(j,pix), im5(j,pix), im6(j,pix), im7(j,pix), im8(j,pix), im9(j,pix)];
% 
%     end
% 
% end
% 
% for k=1:n_sample
% 
%     im0 = double(train0(1200+k,:)');
%     im0 = reshape(im0,28,28);
%     im1 = double(train1(1200+k,:)');
%     im1 = reshape(im1,28,28);
%     im2 = double(train2(1200+k,:)');
%     im2 = reshape(im2,28,28);
%     im3 = double(train3(1200+k,:)');
%     im3 = reshape(im3,28,28);
%     im4 = double(train4(1200+k,:)');
%     im4 = reshape(im4,28,28);
%     im5 = double(train5(1200+k,:)');
%     im5 = reshape(im5,28,28);
%     im6 = double(train6(1200+k,:)');
%     im6 = reshape(im6,28,28);
%     im7 = double(train7(1200+k,:)');
%     im7 = reshape(im7,28,28);
%     im8 = double(train8(1200+k,:)');
%     im8 = reshape(im8,28,28);
%     im9 = double(train9(1200+k,:)');
%     im9 = reshape(im9,28,28);
% 
%     for j=1:28
% 
%         pix = randperm(28,n_pix);
% 
%         I2(j+(k-1)*28,:)=[im0(j,pix), im1(j,pix), im2(j,pix), im3(j,pix), im4(j,pix), im5(j,pix), im6(j,pix), im7(j,pix), im8(j,pix), im9(j,pix)];
% 
%     end
% 
% end
% 
% for k=1:n_sample
% 
%     im0 = double(test0(1+k,:)');
%     im0 = reshape(im0,28,28);
%     im1 = double(test1(1+k,:)');
%     im1 = reshape(im1,28,28);
%     im2 = double(test2(1+k,:)');
%     im2 = reshape(im2,28,28);
%     im3 = double(test3(1+k,:)');
%     im3 = reshape(im3,28,28);
%     im4 = double(test4(1+k,:)');
%     im4 = reshape(im4,28,28);
%     im5 = double(test5(1+k,:)');
%     im5 = reshape(im5,28,28);
%     im6 = double(test6(1+k,:)');
%     im6 = reshape(im6,28,28);
%     im7 = double(test7(1+k,:)');
%     im7 = reshape(im7,28,28);
%     im8 = double(test8(1+k,:)');
%     im8 = reshape(im8,28,28);
%     im9 = double(test9(1+k,:)');
%     im9 = reshape(im9,28,28);
% 
%     for j=1:28
% 
%         pix = randperm(28,n_pix_rhs);
% 
%         %rhs1(j+(k-1)*28,:)=[im0(j,p0), im1(j,p1), im2(j,p2), im3(j,p3), im4(j,p4), im5(j,p5), im6(j,p6), im7(j,p7), im8(j,p8), im9(j,p9)];
%         rhs1(j+(k-1)*28,:)=[im9(j,pix)];
% 
%     end
% 
% end
%
% rhs2 = rhs1;

d = dim_patch;
step = 2;

addpath('./Data')
load mnist_all

n_patch = 13*2;
H1 = zeros(n_patch*dim_patch^2, n_sample*10);
H2 = zeros(n_patch*dim_patch^2, n_sample*10);
I1 = zeros(n_patch*dim_patch^2, n_sample*10);
I2 = zeros(n_patch*dim_patch^2, n_sample*10);


for k=1:n_sample

    im0 = double(train0(1+k,:)');
    im0 = reshape(im0,28,28);
    im1 = double(train1(1+k,:)');
    im1 = reshape(im1,28,28);
    im2 = double(train2(1+k,:)');
    im2 = reshape(im2,28,28);
    im3 = double(train3(1+k,:)');
    im3 = reshape(im3,28,28);
    im4 = double(train4(1+k,:)');
    im4 = reshape(im4,28,28);
    im5 = double(train5(1+k,:)');
    im5 = reshape(im5,28,28);
    im6 = double(train6(1+k,:)');
    im6 = reshape(im6,28,28);
    im7 = double(train7(1+k,:)');
    im7 = reshape(im7,28,28);
    im8 = double(train8(1+k,:)');
    im8 = reshape(im8,28,28);
    im9 = double(train9(1+k,:)');
    im9 = reshape(im9,28,28);

    c = 0;

    for j=1:13
        for i=1:13
            
            c=c+1;

            p0 = im0(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p0 = p0(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2 , k) = p0;
            p1 = im1(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p1 = p1(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample) = p1;
            p2 = im2(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p2 = p2(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*2) = p2;
            p3 = im3(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p3 = p3(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*3) = p3;
            p4 = im4(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p4 = p4(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*4) = p4;
            p5 = im5(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p5 = p5(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*5) = p5;
            p6 = im6(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p6 = p6(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*6) = p6;
            p7 = im7(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p7 = p7(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*7) = p7;
            p8 = im8(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p8 = p8(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*8) = p8;
            p9 = im9(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p9 = p9(:);
            H1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*9) = p9;

        end

    end

end

for k=1:n_sample

    im0 = double(train0(300+k,:)');
    im0 = reshape(im0,28,28);
    im1 = double(train1(300+k,:)');
    im1 = reshape(im1,28,28);
    im2 = double(train2(300+k,:)');
    im2 = reshape(im2,28,28);
    im3 = double(train3(300+k,:)');
    im3 = reshape(im3,28,28);
    im4 = double(train4(300+k,:)');
    im4 = reshape(im4,28,28);
    im5 = double(train5(300+k,:)');
    im5 = reshape(im5,28,28);
    im6 = double(train6(300+k,:)');
    im6 = reshape(im6,28,28);
    im7 = double(train7(300+k,:)');
    im7 = reshape(im7,28,28);
    im8 = double(train8(300+k,:)');
    im8 = reshape(im8,28,28);
    im9 = double(train9(300+k,:)');
    im9 = reshape(im9,28,28);

    c = 0;

    for j=1:13
        for i=1:13
            
            c=c+1;

            p0 = im0(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p0 = p0(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2 , k) = p0;
            p1 = im1(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p1 = p1(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample) = p1;
            p2 = im2(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p2 = p2(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*2) = p2;
            p3 = im3(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p3 = p3(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*3) = p3;
            p4 = im4(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p4 = p4(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*4) = p4;
            p5 = im5(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p5 = p5(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*5) = p5;
            p6 = im6(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p6 = p6(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*6) = p6;
            p7 = im7(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p7 = p7(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*7) = p7;
            p8 = im8(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p8 = p8(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*8) = p8;
            p9 = im9(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p9 = p9(:);
            H2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*9) = p9;

        end
        
    end

end


for k=1:n_sample

    im0 = double(train0(600+k,:)');
    im0 = reshape(im0,28,28);
    im1 = double(train1(600+k,:)');
    im1 = reshape(im1,28,28);
    im2 = double(train2(600+k,:)');
    im2 = reshape(im2,28,28);
    im3 = double(train3(600+k,:)');
    im3 = reshape(im3,28,28);
    im4 = double(train4(600+k,:)');
    im4 = reshape(im4,28,28);
    im5 = double(train5(600+k,:)');
    im5 = reshape(im5,28,28);
    im6 = double(train6(600+k,:)');
    im6 = reshape(im6,28,28);
    im7 = double(train7(600+k,:)');
    im7 = reshape(im7,28,28);
    im8 = double(train8(600+k,:)');
    im8 = reshape(im8,28,28);
    im9 = double(train9(600+k,:)');
    im9 = reshape(im9,28,28);

    c = 0;

    for j=1:13
        for i=1:13
            
            c=c+1;

            p0 = im0(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p0 = p0(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2 , k) = p0;
            p1 = im1(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p1 = p1(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample) = p1;
            p2 = im2(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p2 = p2(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*2) = p2;
            p3 = im3(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p3 = p3(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*3) = p3;
            p4 = im4(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p4 = p4(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*4) = p4;
            p5 = im5(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p5 = p5(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*5) = p5;
            p6 = im6(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p6 = p6(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*6) = p6;
            p7 = im7(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p7 = p7(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*7) = p7;
            p8 = im8(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p8 = p8(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*8) = p8;
            p9 = im9(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p9 = p9(:);
            I1((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*9) = p9;

        end
        
    end

end



for k=1:n_sample

    im0 = double(train0(900+k,:)');
    im0 = reshape(im0,28,28);
    im1 = double(train1(900+k,:)');
    im1 = reshape(im1,28,28);
    im2 = double(train2(900+k,:)');
    im2 = reshape(im2,28,28);
    im3 = double(train3(900+k,:)');
    im3 = reshape(im3,28,28);
    im4 = double(train4(900+k,:)');
    im4 = reshape(im4,28,28);
    im5 = double(train5(900+k,:)');
    im5 = reshape(im5,28,28);
    im6 = double(train6(900+k,:)');
    im6 = reshape(im6,28,28);
    im7 = double(train7(900+k,:)');
    im7 = reshape(im7,28,28);
    im8 = double(train8(900+k,:)');
    im8 = reshape(im8,28,28);
    im9 = double(train9(900+k,:)');
    im9 = reshape(im9,28,28);

    c = 0;

    for j=1:13
        for i=1:13
            
            c=c+1;

            p0 = im0(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p0 = p0(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2 , k) = p0;
            p1 = im1(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p1 = p1(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample) = p1;
            p2 = im2(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p2 = p2(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*2) = p2;
            p3 = im3(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p3 = p3(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*3) = p3;
            p4 = im4(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p4 = p4(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*4) = p4;
            p5 = im5(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p5 = p5(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*5) = p5;
            p6 = im6(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p6 = p6(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*6) = p6;
            p7 = im7(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p7 = p7(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*7) = p7;
            p8 = im8(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p8 = p8(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*8) = p8;
            p9 = im9(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
            p9 = p9(:);
            I2((c-1)*dim_patch^2+1 : c*dim_patch^2, k+n_sample*9) = p9;

        end
        
    end

end



rhs1 = zeros(n_patch*dim_patch^2, 1);

for k=1:n_sample

    t = randi(800);

    im0 = double(train0(t,:)');
    im0 = reshape(im0,28,28);
    im1 = double(train1(t,:)');
    im1 = reshape(im1,28,28);
    im2 = double(train2(t,:)');
    im2 = reshape(im2,28,28);
    im3 = double(train3(t,:)');
    im3 = reshape(im3,28,28);
    im4 = double(train4(t,:)');
    im4 = reshape(im4,28,28);
    im5 = double(train5(t,:)');
    im5 = reshape(im5,28,28);
    im6 = double(train6(t,:)');
    im6 = reshape(im6,28,28);
    im7 = double(train7(t,:)');
    im7 = reshape(im7,28,28);
    im8 = double(train8(t,:)');
    im8 = reshape(im8,28,28);
    im9 = double(train9(t,:)');
    im9 = reshape(im9,28,28);

    c = 0;

    for j=1:13
        for i=1:13
            
            c=c+1;
            
            if digit == 0

                p0 = im0(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p0 = p0(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2 , 1) = p0;

            end

            if digit == 1

                p1 = im1(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p1 = p1(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2, 1) = p1;

            end

            if digit == 2

                p2 = im2(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p2 = p2(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2, 1) = p2;
        
            end

            if digit == 3

                p3 = im3(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p3 = p3(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2, 1) = p3;

            end

            if digit == 4

                p4 = im4(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p4 = p4(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2,1) = p4;

            end

            if digit == 5

                p5 = im5(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p5 = p5(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2, 1) = p5;

            end

            if digit == 6

                p6 = im6(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p6 = p6(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2, 1) = p6;

            end

            if digit == 7

                p7 = im7(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p7 = p7(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2, 1) = p7;

            end

            if digit == 8

                p8 = im8(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p8 = p8(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2, 1) = p8;

            end

            if digit == 9 

                p9 = im9(step*j-d/2+1 : step*j+d/2, step*i-d/2+1 : step*i+d/(2));
                p9 = p9(:);
                rhs1((c-1)*dim_patch^2+1 : c*dim_patch^2, 1) = p9;

            end
            
        end
        
    end

end

rhs2 = rhs1;


end