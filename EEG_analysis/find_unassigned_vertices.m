%%%%written by Kate Christison-Lagay in Feb 2021. Loads the most recent ROI
%%%%map, then finds all vertices that are unassigned, and finds the closest
%%%%ROI and makes them part of it


clear all


toolboxpath = fileparts(which('mni2fs'));
side='lh';
mappeddatadrive='V';
currentfiles=dir (['HOArefinedtable_min200_smoothing6_' side '_' date '*']);

allversions=[];
for n=1:length(currentfiles)
   allfiles=currentfiles(n).name;
   name1=strfind(allfiles, 'version');
   name2=strfind(allfiles, '.mat');
   versionnumber=str2num(allfiles(name1+7:name2-1));
   allversions=[allversions, versionnumber];
   
end

if ~isempty(allversions)
    versionnumber=max(allversions);
    ROIfiles_refined=['HOArefinedtable_min200_smoothing6_' side '_' date '_version' num2str(versionnumber) '.mat'];

else
    versionnumber=0;
    
    previousfiles=dir (['HOArefinedtable_min200_smoothing6_' side '*']);
    filebydate={previousfiles.date}';
    filenamesextracted={previousfiles.name}';
    [extracteddate, orderedfiles]=sortrows(filebydate, 'descend');
    ROIfiles_refined=char(filenamesextracted(orderedfiles(1)));
   
end

   


load allmni2fsfaces.mat

load(ROIfiles_refined)

allfaces_1=['allfaces_' side];
allfaces=eval(allfaces_1);

allvertices_1=['allvertices_' side];
allvertices=eval(allvertices_1);

alldata=refinedtable.Variables;
alldata(:,[1,11])=NaN;

emptyfaces_1=nansum(alldata,2);
emptyfaces=find(emptyfaces_1==0);

for i=1:length(emptyfaces)
    i
    currentface=emptyfaces(i);
    verticesofface=allfaces(currentface,:);
    
    [rindex cindex]=find(allfaces==verticesofface);
    uniquerows=unique(rindex);
    checkassignment1=~ismember(uniquerows, emptyfaces);
    checkassignment=gt(sum(checkassignment1),1);

    initialvertices=allfaces(rindex,:);
    initialvertices2=[initialvertices];
    while ~checkassignment

        all_unique_face_vertices=unique(initialvertices2);
        [rindex cindex]=find(ismember(allfaces, all_unique_face_vertices));
        uniquerows=unique(rindex);
        checkassignment1=~ismember(uniquerows, emptyfaces);
        checkassignment=gt(sum(checkassignment1),1);

        initialvertices=allfaces(rindex,:);
        initialvertices2=[initialvertices2; initialvertices];
        
    end
    
    
    notempties_1=find(checkassignment1==1);
    associatedrows=uniquerows(notempties_1);
    
    nonzerocolumns=find(ge(sum(alldata(associatedrows,:),1),1));
    
    if length(nonzerocolumns>1)
        nonzerocolumn=min(nonzerocolumns);
    else
        nonzerocolumn=nonzerocolumns;
    end
    
    currentarea_name=char(refinedtable.Properties.VariableNames(nonzerocolumn));
    eval(['refinedtable.' currentarea_name '(currentface)=1;']);

end
        

  
    
versionnumber_1=versionnumber+1;
filename=['HOArefinedtable_min200_smoothing6_' side '_' date '_version' num2str(versionnumber_1) '.mat'];
save(filename, 'refinedtable');

    


