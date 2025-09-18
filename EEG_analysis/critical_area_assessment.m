function flag=critical_area_assessment(all_faces, start_vertex, threshold)

%%%function written to determine whether a face belongs to a large enough
%%%area to be included in an ROI. flag of 0 means it passes; flag of 1
%%%means it should be left out of the ROI.

%%written by K Christison-Lagay September 2020


counter=0;
% all_faces=initial_faces;
% % current_face_1=77610;
% current_face_1=all_unique_face_vertices(1);
% threshold=30;
keepingtrack_of_faces=[];
unique_vertices=start_vertex;
unique_vertex=unique_vertices;
currentlyincludedfaces=[];

flag=0;

while lt(counter,threshold)
    
    for a=1:length(unique_vertices)
        
        unique_vertex=unique_vertices(a);
        if ~isempty(find(keepingtrack_of_faces==unique_vertex))
            if isequal(a, length(unique_vertices))
                counter=threshold;
                flag=1;
            end
            
            continue
            
        else
            
            keepingtrack_of_faces=[keepingtrack_of_faces; unique_vertex];
            [rowind colind]=find(all_faces==unique_vertex);
            if isequal(a, length(unique_vertices))&isempty(rowind)
                counter=threshold;
                flag=1;
                continue
            end

            
            currentlyincludedfaces=[currentlyincludedfaces; all_faces(rowind,:)];
            included_faces=unique(currentlyincludedfaces, 'rows');
            counter=length(included_faces);
            
            unique_vertices=unique(included_faces);
            a=1;
            
    
        end
        
        
        
    end
end
