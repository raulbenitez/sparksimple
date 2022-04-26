function [IP] = RegMax3D(IM0,tope)
%pensat per la part positiva del resultat duna cwt normalitzat a 1

steps=51;
%tope=50;% maxnum cands
IM0=double(IM0);

rng=linspace(min(min(min(IM0))),max(max(max(IM0))),steps);

[a0,b0,c0]=size(IM0);
IP=zeros(a0,b0,c0);%

ensenya('3D Localmax. 00%');
contador=0;

for ii=1:steps
    
    thresh=rng(steps-ii+1);
    
    L=bwlabeln(IM0>thresh);
    numcells=max(max(max(L)));
    for jj=1:numcells
        if(numel(unique(IP(L==jj)))==1)% check new object
            [a1,b1]=find(L==jj);
            [b1,c1]=ind2sub([b0,c0],b1);
            aux=round(mean([a1,b1,c1],1));
            contador=contador+1;
            IP(aux(1),aux(2),aux(3))=steps-ii+1;
        end
    end
    if(contador>=tope),break;end
    percentatge(contador,tope+1);
end
percentatge(1,1);




end
