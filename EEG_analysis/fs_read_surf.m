function [surf] = fs_read_surf(fname, v_num, surf_num)
 % fs_read_surf - read a freesurfer surface file
 %
 % [surf] = fs_read_surf(fname)
 %
 % Reads the vertex coordinates (mm) and face lists from a surface file
 % Then finds neighbors for each vertex
 %
 % surf is a structure containg:
 %   nverts: number of vertices
 %   nfaces: number of faces (triangles)
 %   faces:  vertex numbers for each face (3 corners)
 %   coords: x,y,z coords for each vertex
 %   nbrs:   vertex numbers of neighbors for each vertex
 %
 % created:        03/02/06 Don Hagler
 % last modified:  05/09/06 Don Hagler
 %
 % code for reading surfaces taken from Darren Weber's freesurfer_read_surf
 %
 % see also: fs_read_trisurf, fs_find_neighbors, fs_calc_triarea
 %

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 funcname = 'fs_read_surf';

 %QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
 %NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

 TRIANGLE_FILE_MAGIC_NUMBER  =  16777214 ;
 QUAD_FILE_MAGIC_NUMBER      =  16777215 ;

 % open it as a big-endian file
 fid = fopen(fname, 'rb', 'b') ;
 if (fid < 0),
   str = sprintf('%s: could not open surface file %s.',funcname,fname) ;
   error(str) ;
 end

%  magic = fs_fread3(fid) ;
magic=[];
 if (magic == QUAD_FILE_MAGIC_NUMBER),
   surf.nverts = v_num ;
   surf.nfaces = surf_num ;
   fprintf('%s: reading %d quad file vertices...',funcname,surf.nverts);
 tic;
   surf.coords = fread(fid, surf.nverts*3, 'int16') ./ 100 ;
   t=toc; fprintf('done (%0.2f sec)\n',t);
   fprintf('%s: reading %d quad file faces (please wait)...\n',...
     funcname,surf.nfaces); tic;
   surf.faces = zeros(surf.nfaces,4);
   for iface = 1:surf.nfaces,
     for n=1:3
       surf.faces(iface,n) = fs_fread3(fid) ;
     end
     if(~rem(iface, 10000)), fprintf(' %7.0f',iface); end
     if(~rem(iface,100000)), fprintf('\n'); end
   end
   t=toc; fprintf('\ndone (%0.2f sec)\n',t);
 elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER),
   fprintf('%s: reading triangle file...',funcname); tic;
   tline = fgets(fid); % read creation date text line
   tline = fgets(fid); % read info text line

   surf.nverts = fread(fid, 1, 'int32') ; % number of vertices
   surf.nfaces = fread(fid, 1, 'int32') ; % number of faces

   % vertices are read in column format and reshaped below
   surf.coords = fread(fid, surf.nverts*3, 'float32');

   % faces are read in column format and reshaped
   surf.faces = fread(fid, surf.nfaces*3, 'int32') ;
   surf.faces = reshape(surf.faces, 3, surf.nfaces)' ;
   t=toc; fprintf('done (%0.2f sec)\n',t);
 else
   str = sprintf('%s: unknown magic number in surface files.',funcname,fname) ;
   error(str) ;
 end

 surf.coords = reshape(surf.coords, 3, surf.nverts)' ;
 fclose(fid) ;

 %fprintf('...adding 1 to face indices for matlab compatibility.\n\n');
 surf.faces = surf.faces + 1;

 return
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Here is my function for getting neighbors:
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [surf] = fs_find_neighbors(surf)
 % fs_read_surf - find neighboring relations between vertices in a
 surface
 %
 % [surf] = fs_find_neighbors(surf)
 %
 % Input:
 % surf is a structure containg:
 %   nverts: number of vertices
 %   nfaces: number of faces (triangles)
 %   faces:  vertex numbers for each face (3 corners)
 %   coords: x,y,z coords for each vertex
 %
 % Output:
 % surf is a structure containg:
 %   nverts: number of vertices
 %   nfaces: number of faces (triangles)
 %   faces:  vertex numbers for each face (3 corners)
 %   coords: x,y,z coords for each vertex
 %   nbrs:   vertex numbers of neighbors for each vertex
 %
 % created:        05/09/06 Don Hagler
 % last modified:  05/09/06 Don Hagler
 %
 % code for finding neighbors taken from Moo Chung's mni_getmesh
 %
 % see also: fs_read_surf, fs_read_trisurf, fs_calc_triarea
 %

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 funcname = 'fs_find_neighbors';

 if ~isfield(surf, 'faces')
 fprintf('%s: error: input surf must contain faces\n',funcname);
 return;
 end;

 % compute the maximum degree of node -- number of edges = number of
 neighbors
 fprintf('%s: finding number of nearest neighbors...',funcname); tic;
 num_nbrs=zeros(surf.nverts,1);
 for i=1:surf.nfaces
 num_nbrs(surf.faces(i,:))=num_nbrs(surf.faces(i,:))+1;
 end
 max_num_nbrs=max(num_nbrs);
 t=toc; fprintf('done (%0.2f sec)\n',t);

 % find nearest neighbors
 fprintf('%s: finding nearest neighbors...',funcname); tic;
 surf.nbrs=zeros(surf.nverts,max_num_nbrs);
 for i=1:surf.nfaces
 for j=1:3
   vcur = surf.faces(i,j);
   for k=1:3
     if (j ~= k)
       vnbr = surf.faces(i,k);
       if find(surf.nbrs(vcur,:)==vnbr)
         ;
       else
         n_nbr = min(find(surf.nbrs(vcur,:) == 0));
         surf.nbrs(vcur,n_nbr) = vnbr;
       end;
     end;
   end;
 end;
 end;
 t=toc; fprintf('done (%0.2f sec)\n',t);

 return