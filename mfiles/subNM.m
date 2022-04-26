function [h]=subNM(n,m,q,marg)
% subplot fix margin, creates axes in current figure such that it is the
% q-th out of a n-by-m grid with margins specified by marg as a normalised
% value. 
% marg can have two components; first for horizontal second for vertical.
% or four components; left, bot, right, top. 
if(nargin<4),marg=0;end

fila=ceil(q/m);
columna=q-((fila-1)*m);
if(length(marg)==2),marg(3)=marg(1);marg(4)=marg(2);end
if(length(marg)==1),marg(2)=marg(1);marg(4)=marg(1);marg(3)=marg(1);end

%tamany dun sol panell
gs=[(1-(m*marg(1))-marg(3)) (1-(n*marg(2))-marg(4))]./[m n];

pos(1)=(columna*marg(1))+(columna-1)*gs(1);
pos(2)=-gs(2)+1-(((fila-1)*marg(2)+marg(4))+(fila-1)*gs(2));

pos(3)=(columna(end)-columna(1)+1)*gs(1)+(columna(end)-columna(1))*marg(1);
pos(4)=(fila(end)-fila(1)+1)*gs(2)+(fila(end)-fila(1))*marg(2);


h=axes('units','normalized','position',pos);
    




end