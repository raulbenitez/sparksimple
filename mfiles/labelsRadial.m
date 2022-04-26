function [im1]=labelsRadial(pp,im)

[a,b]=size(im);

[my,mx]=find(im>0);
my=round(mean(my));
mx=round(mean(mx));
top=999999;
im0=zeros(size(im));
im1=zeros(size(im));
for ii=1:size(pp,1)


    %vec=pp(ii,2:3)-[my,mx];

    [m,n] = punts2recta(pp(ii,3),pp(ii,2),mx,my);
   %         [top      left   bot   right] 
   coords=zeros(4,2);
    %int top
    x=(1-n)/m;y=1;d=norm([y,x]-pp(ii,2:3));
    if((x<1)||(x>b)),d=top;end
    dists(1)=d;coords(1,:)=[y,x];
    %int bot
    x=(a-n)/m;y=a;d=norm([y,x]-pp(ii,2:3));
    if((x<1)||(x>b)),d=top;end
    dists(3)=d;coords(3,:)=[y,x];
    %int left
    y=m+n;x=1;d=norm([y,x]-pp(ii,2:3));    
    if((y<1)||(y>a)),d=top;end
    dists(2)=d;coords(2,:)=[y,x];
    %int right
    y=m*b+n;x=b;d=norm([y,x]-pp(ii,2:3));
    if((y<1)||(y>a)),d=top;end
    dists(4)=d;coords(4,:)=[y,x];
    
    coords=round(coords);
    
        [aux,qui]=min(dists);
        
     if(aux==top),error('distances shat it');end
    
       T=escriuFrase(num2str(pp(ii,1)),12);T=double(T);T(15:end,:)=[];
     
    [im2] = ellSolet(im0,T,qui,pp(ii,2:3),coords(qui,:),a,b);
    im1=im1+im2;
    if(max(max(im1))>1)
      error('sestan solapant!'); 
    end
end




end


function [im2] = ellSolet(im2,T,qui,p,p2,a,b)
[sv,sh]=size(T);
    mgh=round(sh/2);mgv=round(sv/2);
    
        % creueta
    try im2(p(1)-1:p(1)+1,p(2))=1;end
    try im2(p(1),p(2)-1:p(2)+1)=1;end
    
    switch qui
        case 1;
            x2=p2(2);y2=mgv+2;x2=max([mgh+2,x2]);x2=min([b-sh+mgh,x2]);
        case 2;
            y2=p2(1);x2=mgh+2;y2=max([mgv+2,y2]);y2=min([a-sv+mgv,y2]);
        case 3;
            x2=p2(2);y2=a-sv+mgv-1;x2=max([mgh+2,x2]);x2=min([b-sh+mgh,x2]);
        case 4;
            y2=p2(1);x2=b-sh+mgh-1;y2=max([mgv+2,y2]);y2=min([a-sv+mgv,y2]);
    end
    % numeret
    try im2(y2-mgv:y2-mgv+sv-1,x2-mgh:x2-mgh+sh-1)=T;end
    % 4 marcs
    try im2(y2-mgv:y2-mgv+sv,x2-mgh-1)=1;end
    try im2(y2-mgv:y2-mgv+sv,x2-mgh+sh+1)=1;end
    try im2(y2-mgv-1,x2-mgh-1:x2-mgh+sh+1)=1;end
    try im2(y2-mgv+sv,x2-mgh-1:x2-mgh+sh+1)=1;end
    
    switch qui
        case 1;
            y2=y2-mgv+sv;
        case 2;
            x2=x2-mgh+sh+1;
        case 3;
            y2=y2-mgv-1;
        case 4;
            x2=x2-mgh-1;
    end
    
    % linia fins creueta
    [lesx,lesy]=puntsEnmig([x2,p(2)],[y2,p(1)]);
    for jj=1:length(lesy)
        im2(lesy(jj),lesx(jj))=1;
    end

end