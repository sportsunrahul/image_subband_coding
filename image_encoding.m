% load('hj_optimizer.m');
clc;
clear all;
%% Prepare Filter
% run('hj_optimizer.m')
h0 = [ -0.0323 0.0395  0.0332   -0.0957   -0.0340    0.3321    0.5000    0.3321   -0.0340   -0.0957    0.0332    0.0395   -0.0323];
h1 = [ -0.0323 -0.0395 0.0332    0.0957   -0.0340   -0.3321    0.5000   -0.3321   -0.0340    0.0957    0.0332   -0.0395   -0.0323];
f0 = h0;
f1 = -h1;

%% 1. Read Image
img = imread('lena-gray.bmp');

%% 2. Encoding Image: %%

[L, H] = encode(img, h0, h1);
[LL, LH] = encode(L', h0, h1);
[HL, HH] = encode(H', h0, h1);
LL = LL';
LH = LH';
HL = HL'; 
HH = HH';
LL = uint8(255 * mat2gray(LL));
LH = uint8(255 * mat2gray(LH));
HL = uint8(255 * mat2gray(HL));
HH = uint8(255 * mat2gray(HH));
figure(1);
imshow(img);title("Original Image");
figure(2);
subplot(1,2,1);imshow(uint8(255*mat2gray(L))); title("L");
subplot(1,2,2);imshow(uint8(255*mat2gray(H))); title("H");

figure(3);
sgtitle("Histogram before Quantization");
plot_hist(LL,HL,LH,HH);

%% 3. Quantizing Image
levels = 2^5;
[LL,LH,HL,HH] = quantize(LL,LH,HL,HH,levels);
figure(4);
sgtitle("Histogram after Quantization");
plot_hist(LL,HL,LH,HH);

%% Reconstruction Starts
%% 4. Dequantize Image
[LL,LH,HL,HH] = dequantize(LL,LH,HL,HH,levels);


%% 5. Decode Image
L =  decode(LL', LH', f0, f1);
H =  decode(HL', HH', f0, f1);
L = L';
H = H';
img_new =  decode(L, H, f0, f1);

%% 6. Calulating compression error and displaying final images
img_new = uint8(255 * mat2gray(img_new));
[m,n] = size(img);
l = length(h0);
img_crop = img_new(l:m+l-1,l:n+l-1);
error = ((sum(sum((double(img)-double(img_crop)).^2)))^0.5)/(m*n)
figure(6);
subplot(2,2,1);imshow(img); title("Input Image");
subplot(2,2,2);imshow(img_new); title("Reconstructed Image");
subplot(2,2,3); histogram(img);hold on; histogram(img_new); hold off;
%%

function [L, H] = encode(img, h0, h1)
    [m,n] = size(img);
    l = length(h0);
    TL = zeros(m,n+l-1);
    TH = zeros(m,n+l-1);
    L = zeros(m,ceil(0.5*(n+l-1)));
    H = zeros(m,ceil(0.5*(n+l-1)));
    for i= 1:m
        TL(i,:) = conv(img(i,:),h0);
        TH(i,:) = conv(img(i,:),h1);
        L(i,:) = TL(i,1:2:end);
        H(i,:) = TH(i,1:2:end);
    end  
end
function X = decode(U, V, f0, f1)
    [m,n] = size(U)
    l = length(f0);
    TL = zeros(m,2*n);
    TH = zeros(m,2*n);
    X = zeros(m,2*n+l-1);
    for i = 1:m
        TL(i,1:2:end) = U(i,:); 
        TL(i,2:2:end) = U(i,:);
        TH(i,1:2:end) = V(i,:); 
        TH(i,2:2:end) = V(i,:);
        X(i,:) = conv(TL(i,:),f0);
        X(i,:) = X(i,:) + conv(TH(i,:),f1);
    end
end

function [LL,LH,HL,HH] = quantize(LL,LH,HL,HH,levels)
    b= -62.72;
    a = 9.8e-3;
    LH = a*LH.*LH + b;
    HL = a*HL.*HL + b;
    HH = a*HH.*HH + b;

    LH(LH<0) = 0; LH(LH>255) = 255;
    HL(HL<0) = 0; HL(HL>255) = 255;
    HH(HH<0) = 0; HH(HH>255) = 255;
    figure(5);
    sgtitle("Histogram after Mapping");
    plot_hist(LL,HL,LH,HH);
 
    step_size = 256/levels;
    LL = floor(LL./step_size);
    LH = floor(LH./step_size);
    HL = floor(HL./step_size);
    HH = floor(HH./step_size);
end
function [LL,LH,HL,HH] = dequantize(LL,LH,HL,HH,levels)
    b= -62.72;
    a = 9.8e-3;
    LH = sqrt(double((LH-b)./a));
    HL = sqrt(double((HL-b)./a));
    HH = sqrt(double((HH-b)./a));
    
    step_size = 256/levels;
    LL = LL*step_size;
    LH = LH*step_size;
    HL = HL*step_size;
    HH = HH*step_size;
    
end
function plot_hist(LL,HL,LH,HH);
    subplot(2,2,1);histogram(LL);title('LL');title("LL");
    subplot(2,2,2);histogram(LH);title('LH');title("LH");
    subplot(2,2,3);histogram(HL);title('HL');title("HL");
    subplot(2,2,4);histogram(HH);title('HH');title("HH");
end

