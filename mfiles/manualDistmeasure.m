




% after flos images

for ii=1:length(S)
    res=getres(S(ii).name);
    load([Rfol '/' S(ii).name(1:end-4) '/zData.mat']);
    
I=imread([fol '\' S(ii).name]);
    
     disp([S(ii).name(1:end-4)  ' , ' num2str(estimateZlineDist(pos,res)) 'um']);
    
    
end






for ii=1:length(S)
    res=getres(S(ii).name);
 %   load([Rfol '/' S(ii).name(1:end-4) '/zData.mat']);
    
I=imread([fol '\' S(ii).name]);
figure(555);imshow(I);set(gcf,'position',[ 1          49        2560        1306]);
  [x,y]=ginput(2); 
  
  disp([S(ii).name(1:end-4)  ' , ' num2str(res*norm([diff(x),diff(y)])) 'um']);
%     disp([S(ii).name(1:end-4)  ' , ' num2str(estimateZlineDist(pos,res)) 'um']);
    
    
end

