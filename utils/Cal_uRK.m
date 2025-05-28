function [u0,R0,K0]=Cal_uRK(img0)
u0=sum(img0,2)/(size(img0,2));          
R0=(img0*img0')/(size(img0,2));
K0=((img0-u0)*(img0-u0)')/(size(img0,2));