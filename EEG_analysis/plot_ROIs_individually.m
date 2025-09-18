toolboxpath = fileparts(which('mni2fs'));
side='lh';
mappeddatadrive='V';



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
    ROIfiles=['HOArefinedtable_min200_smoothing6_' side '_' date '_version' num2str(versionnumber) '.mat']

else
    
    previousfiles=dir (['HOArefinedtable_min200_smoothing6_' side '*']);
    filebydate={previousfiles.date}';
    filenamesextracted={previousfiles.name}';
    [extracteddate orderedfiles]=sortrows(filebydate, 'descend');
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


 table_of_faces=refinedtable;
% load blah.mat
for a=1:size(table_of_faces,2)
    figure('Color','w','position',[20 72 800 600])
    
    hold on
    % Load and Render the FreeSurfer surface
    S = [];
    mni2fs_dir = [mappeddatadrive ':\HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master'];
    mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
    mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
    S.hem = side; % choose the hemesphere 'lh' or 'rh'
    S.inflationstep = 5; % 1 no inflation, 6 fully inflated
    S.plotsurf = 'inflated';
    S.lookupsurf = 'mid';
    S.decimation = false; % Decimate the surface for speed. (Use FALSE for publishable quality figures).
    S = mni2fs_brain(S);
    view([-90 0])
    
    
    currentarea=char(table_of_faces.Properties.VariableNames(a))
    currentarea_array=table2array(table_of_faces(:,a));
    set(gcf,'name',currentarea)
    borderexpansions=unique(currentarea_array);
    colortable = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 1 1 1; 0 0 0];

    %%plot 
    for b=2:length(borderexpansions)
        areas_faces=find(currentarea_array(:,1)==(borderexpansions(b)));
        initial_faces=allfaces(areas_faces,:);
        
        S.p=patch('Vertices',S.gfsinf.vertices,'Faces',initial_faces,'EdgeColor','k','EdgeAlpha',0);
        colorlabels = {'y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'};
        cdata = colortable(b-1,:);
        cdata = repmat(cdata,size(initial_faces,1),1);
        Va = ones(size(cdata,1),1).* .8; % can put alpha in here.
        set(S.p,'FaceVertexCData',cdata,'FaceVertexAlphaData',Va,'FaceAlpha',.5);
        shading flat
        
    end
    
%         areas_faces=find(currentarea_array(:,1)~=0);
%         initial_faces=allfaces(areas_faces,:);
%     
%         S.p=patch('Vertices',S.gfsinf.vertices,'Faces',initial_faces,'EdgeColor','k','EdgeAlpha',0);
%         colortable = [1 1 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 1 1 1; 0 0 0];
%         colorlabels = {'y' 'm' 'c' 'r' 'g' 'b' 'w' 'k'};
%         cdata = [1 1 0];
%         cdata = repmat(cdata,size(initial_faces,1),1);
%         Va = ones(size(cdata,1),1).* .8; % can put alpha in here.
%         set(S.p,'FaceVertexCData',cdata,'FaceVertexAlphaData',Va,'FaceAlpha',.5);
%         shading flat
%     
%     
    
end