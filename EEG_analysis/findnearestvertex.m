%%%%Create a rough approximation of ROIs on the other hemisphere given
%%%%refined ROIs of the other hemisphere. Will produce patchy but
%%%%workable ROIs for a jumping off point

%%Written by K Christison-Lagay in March 2021

blah=[];

%%%should load most up to date map:
load HOArefinedtable_min200_smoothing6_rh_19-Feb-2021_version95.mat

load allmni2fsfaces.mat
refinedtablelh=table;
for a=1:size(refinedtable, 2)
    a
    areaname=char(refinedtable.Properties.VariableNames(a));
    currentarea_data=eval(['refinedtable.(''' areaname ''')']);
    currentarea_faces=find(currentarea_data~=0);
    currentvertex=unique(allfaces_rh(currentarea_faces,:));
    allverticesarea=allvertices_rh(currentvertex,:);
    allverticesarea(:,1)=-allverticesarea(:,1);
    left_vertices=ProjectElectrode2TransSurf(size(allvertices_lh,1), allvertices_lh, blah, allverticesarea);
    areafaces=find(ge(sum(ismember(allfaces_lh, left_vertices),2),1));
    
    left_faces=zeros(size(allfaces_lh,1),1);
    left_faces(areafaces,:)=1;

     singletab=table(left_faces, 'VariableNames', {areaname});
     refinedtablelh=[refinedtablelh singletab];
end
versionnumber_1=1;
filename=['HOArefinedtable_min200_smoothing6_lh_' date '_version' num2str(versionnumber_1) '.mat'];
refinedtable=refinedtablelh;
save(filename, 'refinedtable');

function [elecVertices]=ProjectElectrode2TransSurf(vertices_num, vertices, all_electrodes, electrode)
elecVertices = zeros(size(electrode,1),1);
for elec_index = 1:size(electrode,1)
    xCoord = ones(vertices_num,1)*electrode(elec_index,1);
    yCoord = ones(vertices_num,1)*electrode(elec_index,2);
    zCoord = ones(vertices_num,1)*electrode(elec_index,3);
    [minVal, minInd] = min(sqrt((vertices(:,1) - xCoord).^2 + (vertices(:,2) - yCoord).^2 + (vertices(:,3) - zCoord).^2));
    elecVertices(elec_index, 1) = minInd;
end
% load('labels');
% overlay_frame_data = overlay_frame_data';
% electrode_vertex_values = overlay_frame_data(frame_index, all_electrodes);
end