function [I]=pintaFletxa(I,Cul,Cap,col,gruix)

tamany=ceil(size(I,1)*.05)+1;

unitari=(Cap-Cul)/norm(Cap-Cul);

% orientacions cap fletxa
ang=30*3.1416/180;
uP1=[-sin(ang) cos(ang);cos(ang) sin(ang)]*unitari';
ang=-ang;
uP2=[-sin(ang) cos(ang);cos(ang) sin(ang)]*unitari';

% pinta linia
ori=round([Cap]);
des=round(Cul);
    [allx1,ally1] = puntsEnmig([des(2) ori(2)],[des(1) ori(1)]);
    for kk2=1:length(allx1)
        try
        for co=1:3
        I(ally1(kk2):ally1(kk2)+gruix,allx1(kk2):allx1(kk2)+gruix,co)=col(1,co);
        end
        end
    end
%pinta cap    
ori=round(Cap);
des=round([Cap+uP1'*tamany]);
    [allx1,ally1] = puntsEnmig([des(2) ori(2)],[des(1) ori(1)]);
    for kk2=1:length(allx1)
        try
        for co=1:3
        I(ally1(kk2):ally1(kk2)+gruix,allx1(kk2):allx1(kk2)+gruix,co)=col(1,co);
        end
        end
    end
ori=round(Cap);
des=round([Cap+uP2'*tamany]);
    [allx1,ally1] = puntsEnmig([des(2) ori(2)],[des(1) ori(1)]);
    for kk2=1:length(allx1)
        try
        for co=1:3
        I(ally1(kk2):ally1(kk2)+gruix,allx1(kk2):allx1(kk2)+gruix,co)=col(1,co);
        end
        end
    end






end