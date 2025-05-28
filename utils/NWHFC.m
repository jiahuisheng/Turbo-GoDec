function number=NWHFC(HIM,t)
% 
% NWHFC gives the VD number estimated by given false alarm property using
% NWHFC method.
% 
% There are two parameters,NWHFC(HIM,t) where the HIM is the
% Hyperspectral image cube, which is a 3-D data matrix
% [XX,YY,bnd] = size(HIM), XX YY are the image size,
% bnd is the band number of the image cube. 
% t is the false alarm probability.
%
% HFC uses the NWHFC algorithm developed by Dr. Chein-I Chang,
% see http://www.umbc.edu/rssipl/. The Matlab code was
% programmed by Jing Wang in Remote Sensing Signal and 
% Image Processing Lab.
%

[XX,YY,bnd] = size(HIM);
pxl_no = XX*YY;
r = (reshape(HIM,pxl_no,bnd))';

R = (r*r')/pxl_no;
u = mean(r,2);
K = R-u*u';

%======Noise estimation=====
K_Inverse=inv(K); 
tuta=diag(K_Inverse); 
K_noise=1./tuta; 
K_noise=diag(K_noise); 

%=====Noise whitening===
y=inv(sqrtm(K_noise))*r;    
y=reshape(y',XX,YY,bnd);
%=====Call HFC to estimate===
number=HFC(y,t); 
