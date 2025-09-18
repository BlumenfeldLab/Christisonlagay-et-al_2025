function  neighbor_vertex_index=find_neighbor_vertices(faces,vertex_index,step)

% 20190322 Hunki

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
buf_vertex = vertex_index;
neighbor_vertex_index=[];

for step_i=1:step
    buf_neighbor_vertex_index=[];
    for i=1:length(buf_vertex)
        [new_index]=find_neighbor_vertices_from_one_vertex(faces,buf_vertex(i));
        buf_neighbor_vertex_index=[buf_neighbor_vertex_index ; new_index];
    end
    buf_neighbor_vertex_index=unique(buf_neighbor_vertex_index);
    neighbor_vertex_index=[buf_vertex ; buf_neighbor_vertex_index];
    
    buf_vertex=buf_neighbor_vertex_index(~ismember(buf_neighbor_vertex_index,buf_vertex));
end
neighbor_vertex_index=unique(neighbor_vertex_index);

end

function [new_index]=find_neighbor_vertices_from_one_vertex(faces,buf_point)
    
%     [a,~]=find(faces==buf_point);
%     buf_vertices = faces(a,:);
%     new_index=unique(buf_vertices(:));
    
    new_index = unique(faces(faces(:,1) == buf_point | faces(:,2) == buf_point | faces(:,3) == buf_point,:));
    new_index = new_index(new_index~=buf_point);
    
end