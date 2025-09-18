%==================================================================
%
%Title:         CalculateTransparency
%Description:   This function takes the input values from overlay frame
%               data for each electrode and gives it to all the vertices
%               within radius2 distance from each electrode's vertex
%
%Inputs:        "vertices_num"  The number of vertices in the file
%               "elecVertices"  vertex indices of each electrode to be
%                               plotted on this side of the brain
%               "electrode_vertex_values" 
%                               the overlay value assigned to each of those
%                               vertices in elecVertices
%               "overlay"      	[int] this number is an input to the parent
%                               function. If overlay = 0 then it does not
%                               do any weighing based on other electrodes
%                               within a radius2 proximity. All vertices
%                               within radius 2 receive a 1
%
%Outputs:       "vertexCdata"   matrix the length of the number of vertices
%                               with the electrode_vertex_values assigned
%                               to all vertices within radius2 of each
%                               electrode
%==================================================================

%This function should be combined with OverlayColor to speed up processing
%time...

function [ vertexCdata]= CalculateTransparency(vertices_num,vertices, elecVertices, electrode_vertex_values, overlay)
maxAlpha    = 1;
radius1     = 1;
radius2     = 15;
electrodes_num      = length(elecVertices);
electrode_vector    = vertices(elecVertices,:);
vertexCdata = zeros( vertices_num, 1 );
j=1;
k=1;

    for vertIndex = 1 : vertices_num
        sharedTransparency = [];
        sharedElectrode = [];
        
        xCoord = ones(electrodes_num,1)*vertices(vertIndex, 1);
        yCoord = ones(electrodes_num,1)*vertices(vertIndex, 2);
        zCoord = ones(electrodes_num,1)*vertices(vertIndex, 3);

        % Calculate "Transparency" in other words the "weight" based on
        % distance from the vertex
        [distanceToElectrode, electrodeIndice] = sort(sqrt((electrode_vector(:,1) - xCoord).^2 + (electrode_vector(:,2) - yCoord).^2 + (electrode_vector(:,3) - zCoord).^2));

        if overlay == 0
            if any(distanceToElectrode <= radius2)
                vertexCdata(vertIndex,1) = 1;     
            end    
        else

            for n = 1:length(electrodeIndice)
                if distanceToElectrode(n) <= radius1
                    % Calculates the transparencies for each electrode within a
                    % distance less than 15mm from this vertex
                    sharedTransparency = [sharedTransparency; maxAlpha];
                    % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                    sharedElectrode = [sharedElectrode; electrodeIndice(n)]; 
                elseif distanceToElectrode(n) <= radius2
                    % Calculates the transparencies for each electrode within a
                    % distance less than 15mm from this vertex
                    sharedTransparency = [sharedTransparency; (maxAlpha - ((distanceToElectrode(n) - radius1)/(radius2-radius1))*maxAlpha)];
                    % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                    sharedElectrode = [sharedElectrode; electrodeIndice(n)]; 
                end
            end
            % Grabs the zScore values associated with the electrodes that
            % contribute power to this vertex 
            zIn = electrode_vertex_values(sharedElectrode)';
            weightedZ = []; %transparencySum = []; 
            % Calculates the weight of the z score according to the transparencies and zScores associated with this vertex  
            for h = 1:length(zIn)
%                 weightedZ = [weightedZ zIn(h)*(sharedTransparency(h).^2)]
                weightedZ = [weightedZ zIn(h)*(sharedTransparency(h))];
%                 transparencySum = [transparencySum sharedTransparency(h)];
            end
            
            
%             indDelete = []; 
%             for y = 1:length(zIn)
%                 if isnan(weightedZ(y))
%                     indDelete = [indDelete y];
%                 end 
%             end
%             weightedZ(indDelete) = []; 
%             transparencySum(indDelete) = []; 


            weightedZ = nansum(weightedZ); 
%             transparencySum = sum(transparencySum); 
%             weightedZ = weightedZ/transparencySum; 
            % If there is no weighted zScore that means the vertex was greater than
            % 15mm from eevery electrode and has no power or zScore to display in
            % this frame
            if isempty(zIn) 
                vertexCdata(vertIndex,1) = 0;
            elseif isnan(zIn)
                vertexCdata(vertIndex,1) = 0;
            else
                vertexCdata(vertIndex,1) = weightedZ;
            end
        end
    end
end