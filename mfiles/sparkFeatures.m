function [] = sparkFeatures(fol2,Rfol,roiR)

    load([Rfol '/zData0.mat']);
    load([Rfol '/zData1.mat']);
    load([Rfol '/zMetaData.mat']);
    if(exist('mask','var')==1),mk=mask;end
    if(exist('mk','var')==0), mk=cellMask4(sum(volum,3));end
md=bwdist(1-mk);

[DX,DT]=getRes(fol2);
if(DX==0)||(DT==0)
    try
        xmlfil=dir([fol2 '/*.xml']);
        try
            [DX DY DT TTIME NXPIX NYPIX] = physical_parameters(fol2,xmlfil(1)) ;
        catch err
            [DX DY DT TTIME NXPIX NYPIX] = physical_parametersC(fol2,xmlfil(1)) ;
        end
    catch err
        error('Cannot read pixel size and frame rate!');
    end
end




[spkF] = measureFeatures(volum,IX,DX,DT,Rfol,roiR);

if(sum(class(spkF)=='struct')==6)
    
    for ii=1:length(spkF)
        spkF(ii).dist2memb=md(spkF(ii).py,spkF(ii).px)*DX;
    end
else
    spkF= struct([]);
end

%disp(['         Found ' num2str(length(spkF)) ' spark candidates.']);
save([Rfol '/zData2.mat'],'spkF');


end