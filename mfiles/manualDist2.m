function [] = manualDist2()


f=uigetdir;


S=dir([f '/*.xml']);
R=dir([f '/*ch01.tif']);


I=imread([f '/' R.name]);
res=getResolution([f '/' S(1).name])
res=res(1);

figure(555);imshow(I);set(gcf,'position',[ 1          49        2560        1306]);
  [x,y]=ginput(2); 
  
  disp([S(1).name(1:end-4)  ' , ' num2str(res*norm([diff(x),diff(y)])) 'um']);
%     d

end
