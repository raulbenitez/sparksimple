function [f]=deconvTest(g,PSF)





[a1,a2,a3]=size(g);
[b1,b2,b3]=size(PSF);
c=([a1,a2,a3]-[b1,b2,b3])/2;
P2 = padarray(PSF,ceil(c));
dims=ceil(c)==c;
if(dims(1)==0),P2(end,:,:)=[];end
if(dims(2)==0),P2(:,end,:)=[];end
if(dims(3)==0),P2(:,:,end)=[];end
clear PSF;



G=fftn(g);
P=fftn(P2);

%numel(find(P==0))

A=G./P;
f=real((ifftn(A)));




end