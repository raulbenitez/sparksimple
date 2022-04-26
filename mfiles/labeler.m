function [im2] = labeler(im)

[a,b]=size(im);

labs=unique(im);labs(labs==0)=[];

pp=zeros(length(labs),3);
[my,mx]=find(im>0);
my=round(mean(my));
mx=round(mean(mx));

for ii=1:length(labs)

    
    pp(ii,1)=labs(ii);
    [y2,x2]=find(im==labs(ii));
if(~isempty(y2))
   pp(ii,2)=y2(1);pp(ii,3)=x2(1); 
   pp(ii,4)=angol([-1,-1],[x2(1),y2(1)]-[mx,my]);
end
end

try,
    [im2]=labelsRadial(pp,im);
catch  
try
    [im2]=labelsStraight(pp,im);
catch
[im2] = labelsEquidist(pp,im);
end
end

end