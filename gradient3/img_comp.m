close all;
clear all;

% READING THE IMAGE AND CONVERT
IM=imread('lena.jpg');
IM=rgb2gray(IM);
IMD=double(IM);
% SVD DECOMPOSITION OF IMAGE
[U,S,V]=svd(IMD);
% USE DIFFERENT NUMBER OF SINGULAR VALUES TO COMPRESS AND RECONSTRUCT
Derr = [];
num = [];

for N=5:25:200
    % STORE SINGULAR VALUES IN A TEMPORARY VARIABLE
    C = S;
    % DISCARD THE DIAGONAL VALUES NOT REQUIRED FOR COMPRESSION
    C(N+1:end,:)=0;
    C(:,N+1:end)=0;
    % CONSTRUCT IMAGE USING SELECTED SINGULAR VALUES
    D=U*C*V';
    % DISPLAY AND COMPUTE ERROR
    figure;
    buf=sprintf('Image output using %d singular values', N)
    imshow(uint8(D));
    title(buf);
    error=sum(sum((IMD-D).^2));
    % STORE THE VALUES
    Derr = [Derr; error];
    num = [num; N];
end
