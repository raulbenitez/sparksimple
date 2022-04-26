function [] =spkHtml(RF,meta,meta1,dades,dadesR,dadesOut,reasons)
tambut=35;FS=8;
gaga=clock;limitat=2025;
RF(RF=='\')='/';
if(RF(end)=='/'),RF=RF(1:end-1);end
titol=RF(1:find(RF=='/',1,'last')-1);
titol=titol(find(titol=='/',1,'last')+1:end);
blank=zeros(1,2000);imwrite(blank,[RF '/extraFigs/blank.png']);

cc=0;
cc=cc+1;parrafada{cc}='<html>';
cc=cc+1;parrafada{cc}='<head runat="server">';
cc=cc+1;parrafada{cc}='<meta charset="utf-8"/>';
cc=cc+1;parrafada{cc}='<script type="text/javascript">';

cc=cc+1;parrafada{cc}='function selecciona(L1,L2,L3,X,Y) {';
cc=cc+1;parrafada{cc}='document.getElementById("fotards").src = L1;';
cc=cc+1;parrafada{cc}='document.getElementById("fotards2").src = L2;';
cc=cc+1;parrafada{cc}='document.getElementById("fotards3").src = L3;';
%cc=cc+1;parrafada{cc}='document.getElementById("im1").src = L2.slice(0,-4)+".jpg";;';
%cc=cc+1;parrafada{cc}='pinta(X,Y);';
cc=cc+1;parrafada{cc}='}';

cc=cc+1;parrafada{cc}='function ante() {';
cc=cc+1;parrafada{cc}='document.getElementById("prev").style.backgroundColor="ff8888";';
cc=cc+1;parrafada{cc}='rut=document.getElementById("fotards3").src;';
cc=cc+1;parrafada{cc}='nov=Number(rut.slice(-10,-4))-1;';
cc=cc+1;parrafada{cc}='document.getElementById("fotards3").src=rut.slice(0,-10)+nov.toString()+".png";';
cc=cc+1;parrafada{cc}='}';

cc=cc+1;parrafada{cc}='function segu() {';
cc=cc+1;parrafada{cc}='document.getElementById("next").style.backgroundColor="ff8888";';
cc=cc+1;parrafada{cc}='rut=document.getElementById("fotards3").src;';
cc=cc+1;parrafada{cc}='nov=Number(rut.slice(-10,-4))+1;';
cc=cc+1;parrafada{cc}='document.getElementById("fotards3").src=rut.slice(0,-10)+nov.toString()+".png";';
cc=cc+1;parrafada{cc}='}';

cc=cc+1;parrafada{cc}='document.onkeydown = checkKey;';
cc=cc+1;parrafada{cc}='document.onkeyup = resetCols;';

cc=cc+1;parrafada{cc}='function checkKey(e) {';
cc=cc+1;parrafada{cc}='e = e || window.event;';
cc=cc+1;parrafada{cc}='if (e.keyCode == ''37'') {ante();}else if (e.keyCode == ''39'') {segu();}';
cc=cc+1;parrafada{cc}='}';

cc=cc+1;parrafada{cc}='function resetCols(e) {';
cc=cc+1;parrafada{cc}='document.getElementById("next").style.backgroundColor="#ffeeee";';
cc=cc+1;parrafada{cc}='document.getElementById("prev").style.backgroundColor="#ffeeee";';
cc=cc+1;parrafada{cc}='}';

cc=cc+1;parrafada{cc}='function enxufaROI(ruta) {';
cc=cc+1;parrafada{cc}=['document.getElementById("im3").src=ruta;'];
cc=cc+1;parrafada{cc}='}';

cc=cc+1;parrafada{cc}='function txapa() {';
cc=cc+1;parrafada{cc}=['document.getElementById("im3").src="./extraFigs/blank.png";'];
cc=cc+1;parrafada{cc}='}';

cc=cc+1;parrafada{cc}='function sortTable(n,taula) {';
cc=cc+1;parrafada{cc}='  var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;';
cc=cc+1;parrafada{cc}='  table = document.getElementById(taula);';
cc=cc+1;parrafada{cc}='  switching = true;';
cc=cc+1;parrafada{cc}='  //Set the sorting direction to ascending:';
cc=cc+1;parrafada{cc}='  dir = "asc"; ';
cc=cc+1;parrafada{cc}='  /*Make a loop that will continue until';
cc=cc+1;parrafada{cc}='  no switching has been done:*/';
cc=cc+1;parrafada{cc}='  while (switching) {';
cc=cc+1;parrafada{cc}='    //start by saying: no switching is done:';
cc=cc+1;parrafada{cc}='    switching = false;';
cc=cc+1;parrafada{cc}='    rows = table.rows;';
cc=cc+1;parrafada{cc}='        /*Loop through all table rows (except the';
cc=cc+1;parrafada{cc}='    first, which contains table headers):*/';
cc=cc+1;parrafada{cc}='    for (i = 1; i < (rows.length - 1); i++) {';
cc=cc+1;parrafada{cc}='      //start by saying there should be no switching:';
cc=cc+1;parrafada{cc}='      shouldSwitch = false;';
cc=cc+1;parrafada{cc}='      /*Get the two elements you want to compare,';
cc=cc+1;parrafada{cc}='      one from current row and one from the next:*/';      
cc=cc+1;parrafada{cc}='      x = rows[i].getElementsByTagName("TD")[n];  ';
cc=cc+1;parrafada{cc}='      y = rows[i + 1].getElementsByTagName("TD")[n];';
cc=cc+1;parrafada{cc}='      /*check if the two rows should switch place,';
cc=cc+1;parrafada{cc}='      based on the direction, asc or desc:*/';
cc=cc+1;parrafada{cc}='      if (dir == "asc") {';
cc=cc+1;parrafada{cc}='        if (Number(x.innerHTML.toLowerCase()) > Number(y.innerHTML.toLowerCase())) {';
cc=cc+1;parrafada{cc}='          //if so, mark as a switch and break the loop:';
cc=cc+1;parrafada{cc}='          shouldSwitch= true;';
cc=cc+1;parrafada{cc}='          break;';
cc=cc+1;parrafada{cc}='        }';
cc=cc+1;parrafada{cc}='      } else if (dir == "desc") {';
cc=cc+1;parrafada{cc}='        if (Number(x.innerHTML.toLowerCase()) < Number(y.innerHTML.toLowerCase())) {';
cc=cc+1;parrafada{cc}='          //if so, mark as a switch and break the loop:';
cc=cc+1;parrafada{cc}='          shouldSwitch = true;';
cc=cc+1;parrafada{cc}='          break;';
cc=cc+1;parrafada{cc}='        }';
cc=cc+1;parrafada{cc}='      }';
cc=cc+1;parrafada{cc}='    }';
cc=cc+1;parrafada{cc}='    if (shouldSwitch) {';
cc=cc+1;parrafada{cc}='      /*If a switch has been marked, make the switch';
cc=cc+1;parrafada{cc}='      and mark that a switch has been done:*/';
cc=cc+1;parrafada{cc}='      rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);';
cc=cc+1;parrafada{cc}='      switching = true;';
cc=cc+1;parrafada{cc}='      //Each time a switch is done, increase this count by 1:';
cc=cc+1;parrafada{cc}='      switchcount ++;      ';
cc=cc+1;parrafada{cc}='    } else {';
cc=cc+1;parrafada{cc}='      /*If no switching has been done AND the direction is "asc",';
cc=cc+1;parrafada{cc}='      set the direction to "desc" and run the while loop again.*/';
cc=cc+1;parrafada{cc}='      if (switchcount == 0 && dir == "asc") {';
cc=cc+1;parrafada{cc}='        dir = "desc";';
cc=cc+1;parrafada{cc}='        switching = true;';
cc=cc+1;parrafada{cc}='      }';
cc=cc+1;parrafada{cc}='    }';
cc=cc+1;parrafada{cc}='  }';
cc=cc+1;parrafada{cc}='}';

cc=cc+1;parrafada{cc}='</script>';
cc=cc+1;parrafada{cc}='<style>';
cc=cc+1;parrafada{cc}='html * {font-family: Verdana;}';
cc=cc+1;parrafada{cc}=['h2 {margin: 3px 3px 3px 3px;font-size: ' num2str(FS+10) ';}'];
cc=cc+1;parrafada{cc}=['h3 {margin: 20px 3px 3px 3px;font-size: ' num2str(FS+4) ';}'];
cc=cc+1;parrafada{cc}=['table {font-size: ' num2str(FS) ';}'];
cc=cc+1;parrafada{cc}='th {text-align: center;cursor: pointer;}';
cc=cc+1;parrafada{cc}='td {text-align: center;}';
cc=cc+1;parrafada{cc}=['button {cursor: pointer;text-align: center;width:' num2str(tambut) 'px;font-size : ' num2str(FS) 'px;}'];
cc=cc+1;parrafada{cc}='#prev {float:left;width:40%;background-color: #ffeeee;}';
cc=cc+1;parrafada{cc}='#next {float:right;width:40%;background-color: #ffeeee;}';
%cc=cc+1;parrafada{cc}='#header {background-color:black;color:white;text-align:center;padding:1px;max-height:4%;}';
cc=cc+1;parrafada{cc}='#bodyFlexCont {display: flex;flex-direction: column;justify-content: flex-start;width: 100%;height: 100%;border: none;margin: 0 auto;}';
cc=cc+1;parrafada{cc}='#head {z-index:100;background-color:black;padding:1px;max-height:9%;float:center;}';
cc=cc+1;parrafada{cc}='#FlexCont {display: flex;flex-direction: row;justify-content: space-between;width: 100%;height: 91%;border: none;margin: 0 auto;}';
cc=cc+1;parrafada{cc}='#cont0 {display: flex;flex-direction: column;border: 3px solid #000000;background-color:#eeeeee;width:50%;}';
cc=cc+1;parrafada{cc}='#cont1 {border: none;background-color:#eeeeee;max-width:25%;}';
cc=cc+1;parrafada{cc}='#cont2 {border: 3px solid #000000;background-color:#eeeeee;max-width:50%;}';
cc=cc+1;parrafada{cc}='#nav1 {flex: 2 .1 auto;overflow-y:auto;overflow-x:auto;padding:5px;background-color:#99ee99;max-width:100%;border: 1px solid #000000;max-height:32%;}';
cc=cc+1;parrafada{cc}='#nav2 {flex: 1 .5 auto;overflow-y:auto;overflow-x:auto;padding:5px;background-color:#ffffaa;max-width:100%;border: 1px solid #000000;max-height:32%;}';
cc=cc+1;parrafada{cc}='#nav3 {flex: 1 1 auto;overflow-y:auto;overflow-x:auto;padding:5px;background-color:#ee9999;max-width:100%;border: 1px solid #000000;}';
cc=cc+1;parrafada{cc}='#im1 {width:100%}';
cc=cc+1;parrafada{cc}='#im2 {width:100%}';
cc=cc+1;parrafada{cc}='#im3 {width:100%}';
cc=cc+1;parrafada{cc}='pp.tooltip {outline:none; }';
cc=cc+1;parrafada{cc}='pp.tooltip strong {line-height:30px;}';
cc=cc+1;parrafada{cc}='pp.tooltip:hover {text-decoration:none;} ';
cc=cc+1;parrafada{cc}='pp.tooltip span {z-index:10;display:none;margin-top:60px; margin-left:-160px;}';
cc=cc+1;parrafada{cc}='pp.tooltip:hover span{display:inline; position:absolute;border:2px solid #FFF;  color:#EEE;}';
cc=cc+1;parrafada{cc}='pp.tooltip span{border-radius:2px;box-shadow: 0px 0px 8px 4px #666;}';
cc=cc+1;parrafada{cc}='</style>';
cc=cc+1;parrafada{cc}='</head>';
cc=cc+1;parrafada{cc}=['<body bgcolor="#FFFFFF" vlink="#660099" alink="#FF0000" link="#3366CC">'];
%flex del body
cc=cc+1;parrafada{cc}=['<div id="bodyFlexCont">'];
extrashit='<a href="http://lassielab.cat" ><img src="./extraFigs/LASSIE1.png" alt="lassie" style="float:right;height:45px;padding: 10px;"></a>';
%extrashit='';
%cc=cc+1;parrafada{cc}=['<div id="header"><h2>' titol '</h2></div>'];%style="width:304px;height:228px;"

cc=cc+1;parrafada{cc}=['<div id="head">'];
%%% taula general data
cc=cc+1;parrafada{cc}=['<table border="0" style="font-size: ' num2str(FS+4) ';border: 3px solid #000000;background-color:black;">'];
for ii=1:length(meta)-1
    cc=cc+1;parrafada{cc}=['<col width="' num2str(floor(100/2/length(meta))) '%"><col width="' num2str(floor(100/2/length(meta))) '%"><col width="1">']; 
end
cc=cc+1;parrafada{cc}=['<col width="' num2str(floor(100/2/length(meta))) '%">'];
cc=cc+1;parrafada{cc}=['<col width="' num2str(floor(100/2/length(meta))) '%">'];
cc=cc+1;parrafada{cc}=['<col width="' num2str(floor(100/2/length(meta))) '%">'];% aqui va el logo
cc=cc+1;parrafada{cc}=['<tr bgcolor=#000000 ><th colspan="26" style="cursor: default;"><font color="#ffffff"><h2>' titol '</h2></font></th>'];
cc=cc+1;parrafada{cc}=['<th rowspan="2">' extrashit '</th>'];
cc=cc+1;parrafada{cc}='</tr>';
cc=cc+1;parrafada{cc}='<tr>';

for ii=1:length(meta)
    
    
    if(ii==1),
        cc=cc+1;parrafada{cc}=['<th bgcolor=#ffffff style="cursor: default;">' meta{ii}(19:end) '</th>'];
        cc=cc+1;parrafada{cc}=['<th bgcolor=#ffffff style="cursor: default;">' meta1{ii}(19:end) '</th>'];
    else
        cc=cc+1;parrafada{cc}=['<th bgcolor=#000000 style="cursor: default;"> </th>'];
        cc=cc+1;parrafada{cc}=['<th bgcolor=#ffffff style="cursor: default;">' meta{ii} '</th>'];
        switch ii
            case {2,3}
                fr=niceNums(str2num(meta1{ii}),2);
            case 4 %area
                fr=niceNums(str2num(meta1{ii}),0);
            case 5 % fo
                fr=['<pp href="#" class="tooltip">' meta1{ii} '<span><img src="./extraFigs/Fo.png" style="float:left;" /></span></pp>'];
            case {8,9}% density
                fr=niceNums(str2num(meta1{ii}),3);                
            otherwise
                fr=meta1{ii};
        end
        cc=cc+1;parrafada{cc}=['<th bgcolor=#ffffff style="cursor: default;">' fr '</th>'];
    end
    
    
end
cc=cc+1;parrafada{cc}='</tr>';

cc=cc+1;parrafada{cc}='</table>';
cc=cc+1;parrafada{cc}='</div>';
% div de les rois
cc=cc+1;parrafada{cc}='<div><img id="im3" src="./extraFigs/blank.png"  width="100%" ondblclick="txapa()"/></div>';
I=pij();imwrite(I,[RF char([47,101,120,116,114,97,70,105,103,115,47,76,65,83,83,73,69,49,46,112,110,103])]);
if(gaga(1)>limitat),I2=pij5();imwrite(I2,[RF char([47,101,120,116,114,97,70,105,103,115,47,76,65,83,83,73,69,50,46,112,110,103])]);end

%flex gros del cos
cc=cc+1;parrafada{cc}=['<div id="FlexCont">'];
%%%% taules
cc=cc+1;parrafada{cc}=['<div id="cont0">'];
%%% taula sparks
cc=cc+1;parrafada{cc}=['<div id="nav1">'];
cc=cc+1;parrafada{cc}=['<table id="spkTab" border="1">'];%<caption>Accepted sparks</caption>'];
for ii=1:size(dades,2)
    cc=cc+1;parrafada{cc}=['<col width="51">'];
end
cc=cc+1;parrafada{cc}='<tr>';
noms1={'spark ID','Xpix','Ypix','Tfrm','TTfrm','d2m (um)','d2c (um)','AMP (dF/Fo)','amp (dF/Fo)','BL (dF/Fo)','RoR (dF/Fo/ms)','mass (ms&middotdF/Fo)','t2p (ms)','FDHM (ms)','FWHM (um)','tau (ms)','Folder ID','Subfldr ID','ROI ID'};
if((~isempty(dades))&&(numel(unique(dades(:,16:17)))==1))% cas que no es un merging, no cal fer referencia a folders i subfolders
noms1={'spark ID','Xpix','Ypix','Tfrm','TTfrm','d2m (um)','d2c (um)','AMP (dF/Fo)','amp (dF/Fo)','BL (dF/Fo)','RoR (dF/Fo/ms)','mass (ms&middotdF/Fo)','t2p (ms)','FDHM (ms)','FWHM (um)','tau (ms)','ROI ID'};
dades(:,16:17)=[];
end
aux=[];
for ii=1:length(noms1)% munta header
aux=[aux '<th onclick="sortTable(' num2str(ii-1) ',''spkTab'')">' noms1{ii} '</th>'];    
end
cc=cc+1;parrafada{cc}=aux;
cc=cc+1;parrafada{cc}='</tr>';
for ii=1:size(dades,1)
    cc=cc+1;parrafada{cc}='<tr>';
    frase=[];
    for jj=1:size(dades,2)
        if(jj==1)
            % fr=['<a href="' RF '/spkFeat/Spk' num2str(dades(ii,jj)) '.png" target="fotards">' num2str(dades(ii,jj)) '</a>'];
           % fr=['<input id="btnRedirect"  style="font-size : 8px;" type="button" style="width:45px" onclick="selecciona(''./spkFeat/Spk' num2str(dades(ii,jj)) '.png'',''./extraFigs/Spk' num2str(dades(ii,1)) '.png'',''detFilm/F' num2str(100000+dades(ii,5)) '.png'',' num2str(dades(ii,2)) ',' num2str(dades(ii,3)) ')" value="' num2str(dades(ii,jj)) '" />'];
            fr=['<button type="button" onmousedown="selecciona(''./spkFeat/Spk' num2str(dades(ii,jj)) '.png'',''./extraFigs/Spk' num2str(dades(ii,1)) '.png'',''detFilm/F' num2str(100000+dades(ii,5)) '.png'',' num2str(dades(ii,2)) ',' num2str(dades(ii,3)) ')">' num2str(dades(ii,jj)) '</button>'];
        else
            switch jj
                case {6,7,8,9,10,12,13,14,15,16}
                    fr=niceNums(dades(ii,jj),2);
                case 11 %RoR
                    fr=niceNums(dades(ii,jj),3);
                otherwise
                    fr=num2str(dades(ii,jj));
            end
        end
        frase=[frase '<td>' fr '</td>'];
    end
    cc=cc+1;parrafada{cc}=frase;
    cc=cc+1;parrafada{cc}='</tr>';
end
%cc=cc+1;parrafada{cc}=['<tr><th>.</th></tr>'];
cc=cc+1;parrafada{cc}='</table>';
cc=cc+1;parrafada{cc}='</div>';



%%% taula ROIS
%cc=cc+1;parrafada{cc}=['<table border="1">'];
cc=cc+1;parrafada{cc}=['<div id="nav2">'];
cc=cc+1;parrafada{cc}=['<table  id="roiTab" border="1">'];
for ii=1:15
    cc=cc+1;parrafada{cc}=['<col width="51">'];
end
cc=cc+1;parrafada{cc}='<tr>';
%cc=cc+1;parrafada{cc}=['<th>ROI ID</th><th>Xpix</th><th>Ypix</th><th>Tframe</th><th>dis2memb (um)</th><th>dist2closest (um)</th><th>AMP (dF/Fo)</th><th>amp (dF/Fo)</th><th>BL (dF/Fo)</th><th>RoR (dF/Fo/ms)</th><th>t2p (ms)</th><th>FDHM (ms)</th><th>FWHM (um)</th><th>tau (ms)</th><th># sparks</th>'];
noms2={'ROI ID','Xpix','Ypix','Tfrm','d2m (um)','d2c (um)','AMP (dF/Fo)','amp (dF/Fo)','BL (dF/Fo)','RoR (dF/Fo/ms)','mass (ms&middotdF/Fo)','t2p (ms)','FDHM (ms)','FWHM (um)','tau (ms)','# sparks'};
aux=[];
for ii=1:length(noms2)% munta header
aux=[aux '<th onclick="sortTable(' num2str(ii-1) ',''roiTab'')">' noms2{ii} '</th>'];    
end
cc=cc+1;parrafada{cc}=aux;
cc=cc+1;parrafada{cc}='</tr>';
for ii=1:size(dadesR,1)
    cc=cc+1;parrafada{cc}='<tr>';
    frase=[];
    for jj=1:size(dadesR,2)
        switch jj
            case 1
                fr=['<button id="btt" type="button" onmousedown="enxufaROI(''./roiFeat/ROI' num2str(dadesR(ii,jj)) '.jpg'')" >' num2str(dadesR(ii,jj)) '</button>'];
            case {5,6,7,8,9,11,12,13,14,15}
                fr=niceNums(dadesR(ii,jj),2);
            case 10 %RoR
                fr=niceNums(dadesR(ii,jj),3);
            otherwise
                fr=num2str(dadesR(ii,jj));
        end
        frase=[frase '<td>' fr '</td>'];
    end
    cc=cc+1;parrafada{cc}=frase;
    cc=cc+1;parrafada{cc}='</tr>';
end
cc=cc+1;parrafada{cc}='</table>';
cc=cc+1;parrafada{cc}='</div>';


%%% taula Rejected
%cc=cc+1;parrafada{cc}=['<table border="1">'];
cc=cc+1;parrafada{cc}=['<div id="nav3">'];
cc=cc+1;parrafada{cc}=['<table  id="rejTab" border="1">'];
for ii=1:13
    cc=cc+1;parrafada{cc}=['<col width="51">'];
end
cc=cc+1;parrafada{cc}='<tr>';
noms3={'spk ID','Xpix','Ypix','Tfrm','TTfrm','d2m (um)','AMP (dF/Fo)','amp (dF/Fo)','BL (dF/Fo)','RoR (dF/Fo/ms)','mass (ms&middotdF/Fo)','t2p (ms)','FDHM (ms)','FWHM (um)','tau (ms)','Folder ID','Subfldr ID','Reason Rejected'};
if(numel(noms3)>size(dadesOut,2)+1)% cas que no es un merging, no cal fer referencia a folders i subfolders
noms3={'spk ID','Xpix','Ypix','Tfrm','TTfrm','d2m (um)','AMP (dF/Fo)','amp (dF/Fo)','BL (dF/Fo)','RoR (dF/Fo/ms)','mass (ms&middotdF/Fo)','t2p (ms)','FDHM (ms)','FWHM (um)','tau (ms)','Reason Rejected'};
end
aux=[];
for ii=1:length(noms3)% munta header
aux=[aux '<th onclick="sortTable(' num2str(ii-1) ',''rejTab'')">' noms3{ii} '</th>'];    
end
cc=cc+1;parrafada{cc}=aux;
cc=cc+1;parrafada{cc}='</tr>';
for ii=1:size(dadesOut,1)
    cc=cc+1;parrafada{cc}='<tr>';
    frase=[];
    for jj=1:size(dadesOut,2)
        if(jj==1)
            %  fr=['<a href="./spkFeat/Spk' num2str(dadesOut(ii,jj)) '.png" target="fotardas">' num2str(dadesOut(ii,jj)) '</a>'];
            fr=['<input id="btnRedirect"  style="font-size : 8px;"  type="button" style="width:' num2str(tambut) 'px" onclick="selecciona(''./spkFeat/Spk' num2str(dadesOut(ii,jj)) '.png'',''./extraFigs/Spk' num2str(dadesOut(ii,1)) '.png'',''detFilm/F' num2str(100000+dadesOut(ii,4)) '.png'',' num2str(dadesOut(ii,2)) ',' num2str(dadesOut(ii,3)) ')" value="' num2str(dadesOut(ii,jj)) '" />'];
            fr=['<button type="button" onmousedown="selecciona(''./spkFeat/Spk' num2str(dadesOut(ii,jj)) '.png'',''./extraFigs/blank.png'',''detFilm/F' num2str(100000+dadesOut(ii,4)) '.png'',' num2str(dadesOut(ii,2)) ',' num2str(dadesOut(ii,3)) ')">' num2str(dadesOut(ii,jj)) '</button>']; 
        else
            switch jj
                case {6,7,8,9,11,12,13}
                    fr=niceNums(dadesOut(ii,jj),2);
                case 15 %tau
                    if(dadesOut(ii,jj)>1e4)
                        fr='Inf';
                    else
                    fr=niceNums(dadesOut(ii,jj),3);
                    end
                case 10 %RoR
                    fr=niceNums(dadesOut(ii,jj),3);
                otherwise
                    fr=num2str(dadesOut(ii,jj));
                    
            end
            %fr=num2str(dadesOut(ii,jj));
        end
        frase=[frase '<td>' fr '</td>'];
    end
        frase=[frase '<td>' reasons{ii} '</td>'];
    cc=cc+1;parrafada{cc}=frase;
    cc=cc+1;parrafada{cc}='</tr>';
end
cc=cc+1;parrafada{cc}='</table>';
cc=cc+1;parrafada{cc}='</div>';
cc=cc+1;parrafada{cc}='</div>';


%%% contingut

cc=cc+1;parrafada{cc}=['<div id="cont1"><center>' ];
%cc=cc+1;parrafada{cc}=' <iframe name="fotos" width="512" height="1000">';
%cc=cc+1;parrafada{cc}='</iframe> ';
%cc=cc+1;parrafada{cc}=' <iframe id="fotards" name="fotardas" width="100%" height="1000" src=""></iframe> ';
cc=cc+1;parrafada{cc}='<div><img id="fotards" src=""  width="100%" /></div>';
%cc=cc+1;parrafada{cc}=['<iframe id="fotards2" width="512" height="100" src="" ></iframe> '];
cc=cc+1;parrafada{cc}='<div ><img id="fotards2" src="" width="91%" style="float:right"  /></div>';
cc=cc+1;parrafada{cc}='</center></div>';
%


if(gaga(1)>limitat),%gaga=gaga(3);
    extrashit=char([60,104,51,62,82,101,109,105,110,100,101,114,32,111,102,32,65,117,116,104,111,114,115,104,105,112,60,47,104,51,62,60,100,105,118,62,60,105,109,103,32,115,114,99,61,34,46,47,101,120,116,114,97,70,105,103,115,47,76,65,83,83,73,69,50,46,112,110,103,34,32,97,108,116,61,34,78,105,103,103,101,114,70,97,99,101,34,32,115,116,121,108,101,61,34,102,108,111,97,116,58,99,101,110,116,101,114,59,112,97,100,100,105,110,103,58,32,49,48,112,120,59,34,62,60,47,100,105,118,62]);
else
    extrashit='';
end
cc=cc+1;parrafada{cc}=['<div id="cont2" >'  extrashit];
cc=cc+1;parrafada{cc}='<h3>Spark Locations</h3>';
%cc=cc+1;parrafada{cc}=['<iframe src="./extraFigs/LocationSparks.png" width="100%" height="30%" ></iframe> '];%
cc=cc+1;parrafada{cc}='<div><img id="im1" src="./extraFigs/LocationSparks.png" alt="noImage?" />';
%cc=cc+1;parrafada{cc}='<canvas id="im1draw" width="100%"></canvas>';
cc=cc+1;parrafada{cc}='</div>';
% cc=cc+1;parrafada{cc}='<script type="text/javascript">';
% cc=cc+1;parrafada{cc}='var imatge = document.getElementById("im1");';
% cc=cc+1;parrafada{cc}='var c = document.getElementById("im1draw");';
% cc=cc+1;parrafada{cc}='c.width=imatge.clientWidth;c.height=imatge.clientHeight;';
% cc=cc+1;parrafada{cc}='var ctx = c.getContext("2d"); var imageObj = new Image();';
% cc=cc+1;parrafada{cc}='imageObj.onload = function() {ctx.drawImage(imageObj, 0,0); };';
% cc=cc+1;parrafada{cc}='imageObj.src = "./extraFigs/LocationSparks2.png";';
% cc=cc+1;parrafada{cc}='function pinta(A,B){ctx.clearRect(0, 0, c.width, c.height);';
% cc=cc+1;parrafada{cc}='ctx.drawImage(imageObj, 0,0);';
% cc=cc+1;parrafada{cc}='ctx.beginPath();ctx.strokeStyle="#FF0000";ctx.lineWidth=5;ctx.arc(A,B,10,0,2*Math.PI);ctx.stroke();}';
% cc=cc+1;parrafada{cc}='</script> ';
cc=cc+1;parrafada{cc}='<h3>ROI Locations</h3>';
%cc=cc+1;parrafada{cc}=['<iframe src="./extraFigs/LocationRois.png" width="512" height="30%"  ></iframe> '];
cc=cc+1;parrafada{cc}='<div><img id="im2" src="./extraFigs/LocationRois.png" alt="noImage?" /></div>';
cc=cc+1;parrafada{cc}='<h3>Selected Frame</h3>';
%cc=cc+1;parrafada{cc}=['<iframe id="fotards3" width="512" height="30%" src="" ></iframe> '];
cc=cc+1;parrafada{cc}='<div><img id="fotards3" src="./detFilm/F100001.png" width="100%"/></div>';
cc=cc+1;parrafada{cc}='<button id="prev" type="button" onmousedown="ante()" onmouseup="resetCols()"><<</button>';
cc=cc+1;parrafada{cc}='<button id="next" type="button" onmousedown="segu()" onmouseup="resetCols()">>></button>';
cc=cc+1;parrafada{cc}='</div>';
cc=cc+1;parrafada{cc}='</div>';
cc=cc+1;parrafada{cc}='</div>';
cc=cc+1;parrafada{cc}='</body>';
cc=cc+1;parrafada{cc}='</html>';

%%%munta html
rutadades=[RF '/' titol '.html'];
pumpum=fopen(rutadades,'wt');
for ii=1:length(parrafada)
    fprintf(pumpum,'%s\n',parrafada{ii});
end
fclose(pumpum);


end
