function [] = sparkFiltSimple(fol2,Rfol,iT,ta,t2,r2,Bl,Fd)
% no necessita RyR..
ensenya('Filtering Sparks.');
%showDiscardedinfilm=1;
load([Rfol '/zMetaData.mat']);
load([Rfol '/zData0.mat']);
if(yesno==1)
    volum=v1;
    load([Rfol '/zData00.mat']);
    volum=cat(3,volum,v2);
    clear v1;
    clear v2;
end


load([Rfol '/zData2.mat']);
if(exist('ROIs','var')),clear ROIs;end
if(isempty(spkF)==0)
    [DX,DT]=getRes(fol2);
    if(DX==0)||(DT==0),
        xmlfil=dir([fol2 '/*.xml']);
        try
            [DX DY DT TTIME NXPIX NYPIX] = physical_parameters(fol2,xmlfil(1)) ;
        catch err
            [DX DY DT TTIME NXPIX NYPIX] = physical_parametersC(fol2,xmlfil(1)) ;
        end
    end
    filtime=40;%ms de filter window
    w=ceil(filtime/DT);%
    
    IT=sum(volum,3);
    if(exist('mask','var')==0),mask = cellMask4(IT);;end
    %IT=IT.*mask;
    %figure;imagesc(mask)
    IT=(IT-min(min(IT)))/(max(max(IT))-min(min(IT)));
    
    temps=size(volum,3);
    % save([Rfol '/zMetaData.mat'],'fo','mask','temps','DX','DT','roiR');
    
    

    
    [a,b,c]=size(volum);
    
    
    
    
    %MT=double(M>25);gg=MT;
    % gg=double(M>25);gg=bwmorph(gg,'skel',Inf);
    % [MD indist]=bwdist(gg);
    [a0,b0,c0]=size(volum);
    
    Nspks=length(spkF);
    coords=[];dis2ryr=[];
    margeouts=round(0.5/DX);% treure marge!
    margeouts=0;
    dis2ryr=[];
    c2=[];surviving=[];
    
    for ii=1:Nspks
        %  percentatge(ii,Nspks);
        spkF(ii).id=ii;
        spkF(ii).good=0;
        spkF(ii).fail='no';
       % if(aux>w)% inici peli
            if((spkF(ii).px<b0-margeouts+1)&&(spkF(ii).px>margeouts))% per marge
                if((spkF(ii).py<a0-margeouts+1)&&(spkF(ii).py>margeouts))% per marge
                    % algun parametre no mesurable
                    if(sum([isnan(spkF(ii).amp) isnan(spkF(ii).tau) isnan(spkF(ii).ror) isnan(spkF(ii).t2p) isnan(spkF(ii).FDHM) isnan(spkF(ii).AMP) isnan(spkF(ii).BL)])==0)
                        if(spkF(ii).FWHM>0) % FWHM no mesurable
                            if(spkF(ii).BL<=Bl)% baseline
                                if(spkF(ii).amp>=iT)% relative amp
                                    if(spkF(ii).r2>=r2) % decay r2
                                        if(spkF(ii).t2p>=t2(1)) % t2p filter
                                            if(spkF(ii).t2p<=t2(2)) % t2p filter
                                                if(spkF(ii).tau>=ta(1)) % tau filter
                                                    if(spkF(ii).tau<=ta(2)) % tau filter
                                                        if(spkF(ii).FDHM>=Fd(1)) % FDHM filter
                                                            if(spkF(ii).FDHM<=Fd(2)) % FDHM filter
                                                                %thisM=sum(MS(:,:,aux-2:aux),3);
                                                                %surviving=[surviving;ii];
                                                                spkF(ii).good=1;
                                                            else
                                                                spkF(ii).fail='FDHM too large';
                                                            end
                                                        else
                                                            spkF(ii).fail='FDHM too small';
                                                        end
                                                    else
                                                        spkF(ii).fail='tau decay too large';
                                                    end
                                                else
                                                    spkF(ii).fail='tau decay too small';
                                                end
                                            else
                                                spkF(ii).fail='t2p too large';
                                            end
                                        else
                                            spkF(ii).fail='t2p too small';
                                        end
                                    else
                                        spkF(ii).fail='decay r2 too small';
                                    end
                                else
                                    spkF(ii).fail='relative amp too small';
                                end
                            else
                                spkF(ii).fail='previous baseline too large';
                            end
                        else
                            spkF(ii).fail='FWHM unmeasurable';
                        end
                    else
                        spkF(ii).fail='failed parameter (NaN)';
                    end
                else
                    spkF(ii).fail='out of field of view (Y)';
                end
            else
                spkF(ii).fail='out of field of view (X)';
            end
%         else
%             spkF(ii).fail='starting frames';
%         end
    end
    
    
    save([Rfol '/zData2.mat'],'spkF');
    
   
    
    
end

end
