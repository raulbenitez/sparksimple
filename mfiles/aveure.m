function [] = aveure(volum,V3,V4)

figure(88);set(gcf,'position',[86         967        2313         471]);
figure(89);set(gcf,'position',[86         967-471        2313         471]);
figure(90);set(gcf,'position',[86         967-471-471        2313         471]);
temps=size(volum,3);
mx1=max(max(max(volum)));
mx2=max(max(max(V3)));
mx3=max(max(max(V4)));
for ii=1:temps
   figure(88); imagesc(volum(:,:,ii),[0 mx1]); 
   figure(89);imagesc(V3(:,:,ii),[0 mx2]); 
   figure(90);imagesc(V4(:,:,ii),[0 mx3]); 
   pause();
end

end