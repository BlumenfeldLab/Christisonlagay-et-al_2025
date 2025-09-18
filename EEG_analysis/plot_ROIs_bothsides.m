toolboxpath = fileparts(which('mni2fs'));
% close all
    mappeddatadrive='G:\mnt\Data25\';
sides=['rh'; 'lh'];
sidetoplot=[{'rh'}; {'lh'}];



% Load and Render the FreeSurfer surface


for sideinfo=1:length(sides)
    side=sides(sideinfo,:);
    
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
    
    load(ROIfiles)
    
    
    if strcmp(side, 'rh')
        
        refinedtable_rh=refinedtable;
    elseif strcmp(side, 'lh')
        refinedtable_lh=refinedtable;
        
    end
    
end
[combinedareas,ia,ib] = intersect(refinedtable_lh.Properties.VariableNames,refinedtable_rh.Properties.VariableNames, 'stable');
load allmni2fsfaces.mat



allfaces_1=['allfaces_' side];
allfaces=eval(allfaces_1);

plottedareas=listdlg('PromptString', {'Which areas to plot?'}, 'ListString', combinedareas);

% plottedareas=[4 14 36 40 42 43 26];
% plottedareas=1:size(table_of_faces,2);
if isequal(length(plottedareas),size(combinedareas,2))
    colortable = distinguishable_colors(size(combinedareas,2));
else
    %     colortable = distinguishable_colors(size(plottedareas,2));
   % colortable=[0 1 0; 1 .75 0; 1 .48 0; .57 .3 .98;  1 0 1; 0 .86 1;];
    colortable=[1 .75 0];
    % colortable=[.02 .83 .44; .35 .05 .89; .94 .67 .83;  1 .48 0];
    
    % colortable=[0 1 0; 1 .75 0; 0 .4 0; 0 .86 1;  1 .48 0; 1 0 1; .57 .3 .98];
end


for g=1:length(sidetoplot)
    counter=0;
    currentside=sidetoplot(g,:);

    figure('Color','k','position',[20 72 800 600])
    hold on
    % Load and Render the FreeSurfer surface
    S = [];
    mni2fs_dir = [mappeddatadrive 'HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master'];
    mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
    mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
    S.inflationstep = 1; % 1 no inflation, 6 fully inflated, 5 is our default
    S.plotsurf = 'inflated';
    S.lookupsurf = 'inflated';
    S.decimation = false; % Decimate the surface for speed. (Use FALSE for publishable quality figures).
    if strcmp(currentside, 'rh');
        S.hem = 'rh';
        S = mni2fs_brain(S);
        view([90 0])
        
    else
        S.hem = 'lh';
        S = mni2fs_brain(S);
        view([-90 0])
        
    end
    
    for a=1:length(plottedareas) %1:size(table_of_faces,2) %
        counter=counter+1;
        
        
        currentarea=char(combinedareas(plottedareas(a)));
        
        
        
        if strcmp(currentside, 'rh');
            currentarealocation=find(strcmp(refinedtable_rh.Properties.VariableNames, currentarea));
            table_of_faces=refinedtable_rh;
            allfaces=allfaces_rh;
            allvertices=allvertices_rh;
            % choose the hemesphere 'lh' or 'rh'
            
            
            
            
            
        else
            currentarealocation=find(strcmp(refinedtable_lh.Properties.VariableNames, currentarea));
            table_of_faces=refinedtable_lh;
            allfaces=allfaces_lh;
            allvertices=allvertices_lh;
            
        end
        
        
        
        currentarea_array=table2array(table_of_faces(:,currentarealocation));
        
        
        areas_faces=find(currentarea_array(:,1)~=0);
        face_to_vertex1=allfaces(areas_faces,:);
        
        
        
        
        
        
        
        
        
        
        S.p=patch('Vertices',S.gfsinf.vertices,'Faces',face_to_vertex1,'EdgeColor','k','EdgeAlpha',0, 'DisplayName', currentarea);
        cdata =colortable(counter,:);
        cdata = repmat(cdata,size(face_to_vertex1,1),1);
        Va = ones(size(cdata,1),1).* .8; % can put alpha in here.
        set(S.p,'FaceVertexCData',cdata,'FaceVertexAlphaData',Va,'FaceAlpha',.5);
        shading flat
        
        
        
    end
    
    
    legend show % (table_of_faces.Properties.VariableNames)
    hold off;
end