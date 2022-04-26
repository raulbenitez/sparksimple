function [C] = conv3(A,B,verbose)
% Computes 3D convolution using FFT.
% Returns the same as if using convn(A,B,'same') and assumes size(A)>=size(B).
% Argument verbose is optional to display progress.


if(nargin<2),error('Expecting 2 arguments');end
if(nargin==2),verbose=0;end
        
if(verbose==1),fprintf(1,'\b Resizing volumes.');end
[a1,a2,a3]=size(A);
[b1,b2,b3]=size(B);
c=([a1,a2,a3]-[b1,b2,b3])/2;
B2 = padarray(B,ceil(c));
dims=ceil(c)==c;
if(dims(1)==0),B2(end,:,:)=[];end
if(dims(2)==0),B2(:,end,:)=[];end
if(dims(3)==0),B2(:,:,end)=[];end
clear B;

if(verbose==1),fprintf(1,' FFT 1.');end
AA=fftn(A);clear A;

if(verbose==1),fprintf(1,' FFT 2.');end
BB=fftn(B2);clear B2;

if(verbose==1),fprintf(1,' Product.');end
CC=AA.*BB;
clear AA;
clear BB;
if(verbose==1),fprintf(1,' iFFT.');end
C=real(fftshift(ifftn(CC)));
clear CC;

% if(verbose==1),fprintf(1,' Resizing result.\n');end
% C = padarray(C,[1 1 1]);
% C(end-1:end,:,:)=[];C(:,end-1:end,:)=[];C(:,:,end-1:end)=[];



end