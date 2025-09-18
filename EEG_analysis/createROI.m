%%%%written by Kate Christison-Lagay in fall 2020. Given all the faces and
%%%%vertices for a brain mesh, and given a specific set of coordinates, it
%%%%creates an ROI of a given size using a nearest neighbors method


function outputarray=createROI(allfaces, allvertices, xcoor, ycoor, zcoor, targetsize)
xes=(lt(abs(allvertices(:,1)-xcoor),.01));
yes=(lt(abs(allvertices(:,2)-ycoor),.01));
zes=(lt(abs(allvertices(:,3)-zcoor),.01));
seedvertex=find(sum([xes, yes, zes], 2)==3);

% blah=[];
% 
% seedvertex=ProjectElectrode2TransSurf(size(allvertices,1), allvertices, blah, [xcoor, ycoor, zcoor]);


current_unique_face_vertex=seedvertex;
%%now go through, find all vertices, and if its neighbors aren't
%%added, add them.
outputarray=zeros(length(allfaces),1);
[rindex cindex]=find(allfaces==current_unique_face_vertex);

initial_faces=allfaces(rindex,:);
while lt(size(initial_faces,1), targetsize)
    all_unique_face_vertices=unique(initial_faces);
        for a=1:length(all_unique_face_vertices)
            current_unique_face_vertex=all_unique_face_vertices(a);
            %%now go through, find all vertices, and if its neighbors aren't
            %%added, add them.
            [rindex cindex]=find(allfaces==current_unique_face_vertex);
            if gt(length(rindex),4)
                to_add=allfaces(rindex,:);
                initial_faces=[initial_faces; to_add];
                outputarray(rindex)=1;
            end
        end
        
end
%     outputarray(rindex)=outputarray(rindex)+((output_array(rindex)==0)*(roisnumber+(smooththis*.1)));
end




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


% function outputvertices=createROI(seedvertex, size)
%
% current_unique_face_vertex=seedvertex;
% %%now go through, find all vertices, and if its neighbors aren't
% %%added, add them.
% [rindex cindex]=find(S.gfsinf.faces==current_unique_face_vertex);
% if gt(length(rindex),4)
%     to_add=S.gfsinf.faces(rindex,:);
%     initial_faces=[initial_faces; to_add];
%     output_array(rindex)=output_array(rindex)+((output_array(rindex)==0)*(roisnumber+(smooththis*.1)));
% end
%
%
%
% end