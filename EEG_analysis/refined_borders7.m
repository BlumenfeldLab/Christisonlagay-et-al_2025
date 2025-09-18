%%%%written by Kate Christison-Lagay in Feb 2021. Loads the most recent ROI
%%%%map, and you can specify which areas should be added together to form a
%%%%combined area (and gives this area a new name). Expects to be
%%%%fed an ROI file that has 'refinedtable' as one of its variables, but it
%%%%dynamically updates the ROIs that are presented.


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
    filebydate2=datetime(filebydate, 'InputFormat', 'dd-MMMM-yyyy HH:mm:ss');
    filenamesextracted={previousfiles.name}';
    
    [extracteddate orderedfiles]=sortrows(filebydate2, 'descend');
    ROIfiles_refined=char(filenamesextracted(orderedfiles(1)));
    
end




load allmni2fsfaces.mat

load(ROIfiles_refined)

allfaces_1=['allfaces_' side];
allfaces=eval(allfaces_1);

allvertices_1=['allvertices_' side];
allvertices=eval(allvertices_1);


% table_converted=table2array(table_of_facesHOA);
% table_converted_1=table_converted;

%%%

ROIs_1={refinedtable.Properties.VariableNames};
ROIareas=ROIs_1{1}';
f=listdlg('PromptString', {'What areas to add together?'}, 'ListString', ROIareas);

if (length(f)>1)
    newareaname1=inputdlg('What should the new combined area be named?');
    newareaname=strrep(newareaname1, ' ', '');
    while ~isempty(find(strcmp(ROIareas, newareaname)))
        
        newareaname1=inputdlg('That name already exists. Please enter new name:');
        newareaname=strrep(newareaname1, ' ', '');
    end
    
    
    faces_combined_1=[];
    
    for i=1:length(f)
        currentselection=f(i);
        currentarea_name=char(ROIareas(currentselection));
        
        eval(['existingarea=find(refinedtable.' currentarea_name '~=0);'])
        
        faces_combined_1=[faces_combined_1; existingarea];
        
        
    end
    
    finalfaces=unique(faces_combined_1); 
    faces_area=zeros(size(allfaces,1),1);
    faces_area(finalfaces,:)=1;
    newareaname2=char(newareaname);
    eval(['refinedtable.' newareaname2 '=faces_area;']);
    
    versionnumber_1=versionnumber+1;
    filename=['HOArefinedtable_min200_smoothing6_' side '_' date '_version' num2str(versionnumber_1) '.mat'];
    save(filename, 'refinedtable');
    
    
    
    
    
end

