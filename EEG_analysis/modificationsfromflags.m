%%%%written by Kate Christison-Lagay in Feb 2021. requires a brain to be
%%%%open that has flagged datapoints; then gets fed how large of a modification
%%%%size (using a nearest neighbors method), and should be given allfaces
%%%%and allvertices from the current brain mesh (at least for KCL:mni2fs's
%%%%version).

function [finalfaces] = modificationsfromflags(flaggedpoints, modificationsize, allfaces, allvertices)

flaggedpoints1=flaggedpoints{1};

if isequal(length(modificationsize),1)
    modificationsizes=modificationsize.*ones(size(flaggedpoints1, 1),1);
elseif isequal(length(modificationsize) ,size(flaggedpoints1, 1))
    
    modificationsizes=modificationsize;
else
    error('Mismatch in number of flagged points and modification size!')
end
    
faces_combined_1=[];
for i=1:size(flaggedpoints1, 1)
    modsize=modificationsizes(i);
    xcoordinate_1=flaggedpoints1{i}{1};
    ycoordinate_1=flaggedpoints1{i}{2};
    zcoordinate_1=flaggedpoints1{i}{3};
    
    xcoordinate=extractNumFromStr(xcoordinate_1);
    ycoordinate=extractNumFromStr(ycoordinate_1);
    zcoordinate=extractNumFromStr(zcoordinate_1);
    
    
    faces1_1=createROI(allfaces, allvertices, xcoordinate, ycoordinate, zcoordinate, modsize);
%     faces1=find(faces1_1~=0);
    faces_combined_1=[faces_combined_1, faces1_1];
end
finalfaces=find(sum(faces_combined_1,2)~=0);


    function numArray = extractNumFromStr(str)
        str1 = regexprep(str,'[,;=]', ' ');
        str2 = regexprep(regexprep(str1,'[^- 0-9.eE(,)/]',''), ' \D* ',' ');
        str3 = regexprep(str2, {'\.\s','\E\s','\e\s','\s\E','\s\e'},' ');
        numArray = str2num(str3);
    end



end