%%%Find the faces of an ROI, given the vertices of that ROI
%%% Written by K Christison-Lagay in March 2021

function [tableoutput]=findfacesfromvertices(faces, verticesofinterest, nameofarea)
[rindex cindex]=find(ismember(faces, verticesofinterest));

faces_area=zeros(size(faces,1),1);
faces_area(rindex,:)=1;

tableoutput=table(faces_area, 'VariableNames', { nameofarea });
% save(filename, 'refinedtable');
end