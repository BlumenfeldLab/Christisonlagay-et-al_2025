toolboxpath = fileparts(which('mni2fs'));
% close all
mappeddatadrive='V';

side='rh';


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
    ROIfiles=['HOArefinedtable_min200_smoothing6_' side '_' date '_version' num2str(versionnumber) '.mat']

else
    
    previousfiles=dir (['HOArefinedtable_min200_smoothing6_' side '*']);
    filebydate={previousfiles.date}';
    filebydate2=datetime(filebydate, 'InputFormat', 'dd-MMMM-yyyy HH:mm:ss');
    filenamesextracted={previousfiles.name}';
    
    [extracteddate orderedfiles]=sortrows(filebydate2, 'descend');
    ROIfiles=char(filenamesextracted(orderedfiles(1)));
   
end




load allmni2fsfaces.mat

load(ROIfiles)
allfaces_1=['allfaces_' side];
allfaces=eval(allfaces_1);

if contains(ROIfiles, 'HOAROI')
    table_of_faces=table_of_facesHOA;
elseif contains(ROIfiles, 'refined')
    table_of_faces=refinedtable;
end

figure('Color','k','position',[20 72 800 600])
hold on
% Load and Render the FreeSurfer surface
S = [];
mni2fs_dir = [mappeddatadrive ':\HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master'];
mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
S.hem = side; % choose the hemesphere 'lh' or 'rh'
S.inflationstep = 5; % 1 no inflation, 6 fully inflated, 5 is our default
S.plotsurf = 'inflated';
S.lookupsurf = 'mid';
S.decimation = false; % Decimate the surface for speed. (Use FALSE for publishable quality figures).
S = mni2fs_brain(S);
% mni2fs_lights('on', [.8 .2 .1])
view([-90 0])

ROIs_1={table_of_faces.Properties.VariableNames};
ROIareas=ROIs_1{1}';
plottedareas=listdlg('PromptString', {'What areas to add together?'}, 'ListString', ROIareas);


% plottedareas=[4 14 36 40 42 43 26];
% plottedareas=1:size(table_of_faces,2);
if isequal(length(plottedareas),size(table_of_faces,2))
    colortable = distinguishable_colors(size(table_of_faces,2));
else
%     colortable = distinguishable_colors(size(plottedareas,2));
colortable=[0 1 0; 1 .75 0; 1 .48 0; .57 .3 .98;  1 0 1; 0 .86 1;];
% colortable=[.02 .83 .44; .35 .05 .89; .94 .67 .83;  1 .48 0];

% colortable=[0 1 0; 1 .75 0; 0 .4 0; 0 .86 1;  1 .48 0; 1 0 1; .57 .3 .98];
end
            counter=0;

    for a=plottedareas %1:size(table_of_faces,2) % 
            counter=counter+1;

      %  [18, 20, 21, 33, 36, 30, 4, 19, 12, 8, 10, 16, 17, 11, 14, 15, 34, 32, 1, 2, 3] %
        %     areas_faces=find(outputarray~=0);
        % areas_faces=supfront_faces;
        currentarea=char(table_of_faces.Properties.VariableNames(a))
        currentarea_array=table2array(table_of_faces(:,a));
        
        %     borderexpansions=unique(currentarea_array);
        
        %
        areas_faces=find(currentarea_array(:,1)~=0);
        initial_faces=allfaces(areas_faces,:);
        
        S.p=patch('Vertices',S.gfsinf.vertices,'Faces',initial_faces,'EdgeColor','k','EdgeAlpha',0, 'DisplayName', currentarea);
        cdata =colortable(counter,:);
        cdata = repmat(cdata,size(initial_faces,1),1);
        Va = ones(size(cdata,1),1).* .8; % can put alpha in here.
        set(S.p,'FaceVertexCData',cdata,'FaceVertexAlphaData',Va,'FaceAlpha',.5);
        shading flat
        
        %     legend show
        
        
    end


 legend show % (table_of_faces.Properties.VariableNames)