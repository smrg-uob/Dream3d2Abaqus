function dream2abq(voxFileName, inpFileName)
    %creates an input file containing periodic boundary conditions using
    %the *dream2abhex* conversion approach developed by:
        %Marat I. Latypov while at POSTECH
        %latmarat@postech.edu
        %Read more at http://latmarat.net/blog/scripts/dream2abahex/
    %the Dream.3D conversion reads the file generated using the exporting
    %to Los Alamos FFT filter (extension: .vox) and generates:
        % - element sets with individual phases
        % - sections with individual phases
        % - element sets with grains the nodes
        % - euler angles for each grain
    %The Abaqus input file creation has been updated to also now include:
        % - the construction of material files for each grain using the
            %orientation of the grain and the rotation matrix
        % - periodic boundary conditions which has currently been developed
            % for an axial loading in the x-direction
    
    %% The Syntax used to create the input:
        %vox2aba(voxFileName, inpFileName)
  
    %% Input
        % voxFileName - full path to Los Alamos FFT file written by Dream.3D
        % inpFileName - desired path to Abaqus input file 
            % Example:
            % vox2aba('dp_64x64x64.vox', 'dp_64x64x64-aba.inp')
        %material parameter file which contains one sheet for the desired 
        %parameters to be assigned to the grains:
            %Sheet Name: Material_parameters
        %centroid location (x,y,z):
            %Sheet name: centroid
        %The name of the excel file which contains these sheets should be
        %named the following:
            %inputfile_info.xlsx
    % % Result
    % ABAQUS input file containing C2D8 mesh with 
    % - element sets with individual phases
    % - sections with individual phases
    % - element sets with grains 
    % - node sets of faces for easier assignment of BCs
    % and orientations (tex) file containing
    % Euler angle sets for each element in the mesh
    %
    % Read more at http://latmarat.net/blog/scripts/dream2abahex/
   
    % --------------------------
    % written by
    %Dylan Agius July 2018
    %dylan.agius@bristol.ac.uk
    
    % --------------------------
    
    tic  
    %% Digest data from Los Alamos FFT file generated in Dream.3D
    format shortg
    % open FFT file
    fid = fopen(voxFileName,'rt');
    rawData = textscan(fid, '%f %f %f %f %f %f %d %d','delimiter',' ');
    fclose(fid);

    % load euler angles, coordinates, grain and phase IDs
    euler   = cell2mat(rawData(1:3));
    xyz     = cell2mat(rawData(4:6));
    grains  = cell2mat(rawData(7));
    phases  = cell2mat(rawData(8));
    [~,order] = sortrows(xyz,[1,3,2]);
    
    %grains=grains+1;
    
    clear max_value
    clear total_els
    %calculate the maximum value of elements and nodes
    max_value=max(cell2mat(rawData(4)));
    total_els=max_value^3;
    total_nodes=(max_value+1)^3;
    
    %coupling node
    middle_nodex=(max_value/2);
    middle_nodey=0.;
    middle_nodez=(max_value/2);
    
    %grains and corresponding euler angles
    grain_record(1)=grains(1);
    g=1;
    sort(1)=0;
    ratio(1)=0;
    for i=1:length(grains)
       %%
       clear ratio
       if(i+1<=length(grains))
           %%
           grain_record(i)=grains(i);
           if(grain_record(i)~=grains(i+1))
                %%
                for l=1:length(sort)
                    ratio(l)=(double(sort(l))/double(grain_record(i)));
                end
                %%
                if(ismember([1],double(ratio))==0)
                    euler_angle1(g)=euler(i,1);
                    euler_angle2(g)=euler(i,2);
                    euler_angle3(g)=euler(i,3);
                 
                    grain_order(g)=grain_record(i);
                    phase_order(g)=phases(i);
                    
                    g=g+1;
                    sort(g)=grain_record(i);
                end 
           end
       end
    end
   

   
   
    %% Generate the rotation matrix        
    for ii=1:length(grain_order)
        %%
        %need to create the rotation matrix for each euler angle for each
        %grain.  This matrix is then used to rotate a global orientation to
        %the what is is n the local orientation.
        zrot=[cosd(euler_angle1(ii)), sind(euler_angle1(ii)), 0; -sind(euler_angle1(ii)), cosd(euler_angle1(ii)),0; 0,0,1];
        xrot=[1,0,0;0,cosd(euler_angle2(ii)),sind(euler_angle2(ii));0,-sind(euler_angle2(ii)),cosd(euler_angle2(ii))];
        zrot2=[cosd(euler_angle3(ii)),sind(euler_angle3(ii)),0;-sind(euler_angle3(ii)),cosd(euler_angle3(ii)),0;0,0,1];
   
        %total rotation matrix
        total_rot=transpose(zrot2*xrot*zrot);
   
        %rotating the different vectors from the local grain system to the 
        %global system.
        %first considering the vector [0,0,1]
        v1=[0;0;1];
        rot_v1=total_rot*v1;
   
        u(grain_order(ii))=rot_v1(1,1);
        v(grain_order(ii))=rot_v1(2,1);
        w(grain_order(ii))=rot_v1(3,1);
   
        %next considering the vector [0,1,0]
        v2=[0;1;0];
        rot_v2=total_rot*v2;
   
        h(grain_order(ii))=rot_v2(1,1);
        k(grain_order(ii))=rot_v2(2,1);
        l(grain_order(ii))=rot_v2(3,1);
  

    end
         
    % get the number of vox along x, y, z
    xVox = size(unique(xyz(:,1)),1);
    yVox = size(unique(xyz(:,2)),1);
    zVox = size(unique(xyz(:,3)),1);

    % get step size and boundaries for the mesh
    step   = zeros(1,3);
    boxmin = zeros(1,3);
    boxmax = zeros(1,3);
    for ii = 1:3
        %%
        step(ii) = min(diff(unique(xyz(:,ii))));
        boxmin(ii) = xyz(1,ii)-step(ii)/2;
        boxmax(ii) = xyz(end,ii)+step(ii)/2;
    end

    %% Generate 3D mesh
    % generate nodes 
    [x,y,z] = meshgrid(boxmin(1):step(1):boxmax(1),boxmin(2):step(2):boxmax(2),boxmin(3):step(3):boxmax(3));
    numNodes = numel(x);
    coord = [reshape(x,numNodes,1), reshape(y,numNodes,1), reshape(z,numNodes,1)];
    nodes = [(1:numNodes)', sortrows(coord,[1,3,2])];

    % allocate array for elements
    elem = zeros(size(xyz,1),9);
    count = 1;

    % start loop over voxel dimensions
    for ix = 1:xVox
        %%
        for iz = 1:zVox
            %%
            for iy = 1:yVox
                %%
                % get element label
                elem(count,1) = count;

                % nodes on the plane with lower x
                elem(count,2) = iy + (iz-1)*(yVox+1) + (ix-1)*(yVox+1)*(zVox+1);
                elem(count,3) = elem(count,2) + 1;
                elem(count,4) = elem(count,3) + yVox + 1;
                elem(count,5) = elem(count,2) + yVox + 1;

                % nodes on the plane with higher x
                elem(count,6) = iy + (iz-1)*(yVox+1) + ix*(yVox+1)*(zVox+1);
                elem(count,7) = elem(count,6) + 1;
                elem(count,8) = elem(count,7) + yVox + 1;
                elem(count,9) = elem(count,6) + yVox + 1;

                count = count+1;
            end
        end
    end

    %% Write the overall element and node sets to input file
    % open inp file and write keywords 
    inpFile = fopen(inpFileName,'wt');
    fprintf(inpFile,'** Generated by: dream2abahex.m\n');
    fprintf(inpFile,'**PARTS\n**\n');
    fprintf(inpFile,'*Part, name=DREAM\n');

    % write nodes
    fprintf(inpFile,'*NODE\n');
    fprintf(inpFile,'%d,\t%e,\t%e, \t%e\n',nodes');

    % write elements
    fprintf(inpFile,'*Element, type=C3D8\n');
    fprintf(inpFile,'%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d\n',elem');
    %% Write the elements sets for each grain to the input file
    % create element sets containing grains
    for ii = 1:numel(unique(grains))
        %%
        fprintf(inpFile,'\n*Elset, elset=GRAIN-%d\n',grain_order(ii));
        fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',elem(grains==grain_order(ii))');
        numels=0;
        
        for tt=1:length(elem(grains==grain_order(ii)))
            %%
            numels=numels+1;
        end
       numels_total(grain_order(ii))=numels;
    end
    %% Write element set for each phase to input file
    uniPhases = unique(phases);
    for ii = 1:numel(unique(phases))
        fprintf(inpFile,'\n*Elset, elset=Phase-%d\n',ii);
        fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',elem(phases==uniPhases(ii))');
    end
    %% Calculate grain spherical equivalent diameter
    % calulate diamater in microns
    % additionally, the dimaters for each ground are written to a separate
    % text file to be used to developed a grain size histogram
    diameterID=fopen('diameter.txt','w');
    for ii=1:numel(unique(grains))
        %%
        diameter(grain_order(ii))=((((6.0/pi)*(numels_total(grain_order(ii))))^(1/3)));
        fprintf(diameterID, '%d\n', diameter(grain_order(ii)));
    end
    fclose(diameterID);
        
    %% write sections to each grain
    for ii=1:length(grain_order)
        %%
        fprintf(inpFile,'\n**Section: Section_Grain-%d\n*Solid Section, elset=GRAIN-%d, material=MATERIAL-GRAIN%d\n,\n',grain_order(ii),grain_order(ii),grain_order(ii));
    end
    %% Continue writing the input file with assembly information
    % write a closing keyword
    fprintf(inpFile,'*End Part');
    
    %writing assembly
    fprintf(inpFile,'\n**\n**ASSEMBLY\n**');
    fprintf(inpFile,'\n*Assembly, name=Assembly\n**');
    fprintf(inpFile,'\n*Instance, name=DREAM-1, part=DREAM\n');
    
    % write nodes
    fprintf(inpFile,'*NODE\n');
    fprintf(inpFile,'%d,\t%e,\t%e, \t%e\n',nodes');

    % write elements
    fprintf(inpFile,'*Element, type=C3D8\n');
    fprintf(inpFile,'%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d,%8d\n',elem');
    
    fprintf(inpFile,'\n*End Instance\n**');
    
    %% Create the reference nodes the displacement/force controls are applied.
    fprintf(inpFile,'\n*Node\n1,%f,%f,%f\n',middle_nodex,0,-1);
    fprintf(inpFile,'\n*Node\n2,%f,%f,%f\n',-1,0,middle_nodez);
    fprintf(inpFile,'\n*Node\n3,%f,%f,%f\n',middle_nodez,-1,middle_nodez);
    fprintf(inpFile,'*Nset, nset=refnodez\n 1,\n');
    fprintf(inpFile,'*Nset, nset=refnodex\n 2,\n');
    fprintf(inpFile,'*Nset, nset=refnodey\n 3,\n');
    
   %% Applying periodic boundary conditions
   % Node sets defining the faces, edges and veritices of the RVE are
   % developed in the following for loops.  The sets defining a face do not
   % include the nodes at the edges because they are defined in a separate
   % node set.  Additionally, the node sets defining edges did not contain
   % vertices since they were their own seperate node set.
  
   %left face z-axis
   kn=1;
   total=1;
   total_all=1;
   for jn=1:(max_value+1)-2
       %%
       for in=kn:(((max_value+1)-3)+kn)
           %%
           set_node_rights(total)=(((max_value+1)^2)+1)+in;
           fprintf(inpFile,'\n*Nset, nset=RightS-%d, instance=DREAM-1\n',total);
           fprintf(inpFile,'%d,\n',set_node_rights(total)');
           total=total+1;
       end
       total_all=total_all+total;
       kn=kn+(max_value+1)^2;
   end
  
   %right face z-axis
   kn=((max_value+1)*(max_value));
   total=1;
   total_all=1;
   for jn=1:(max_value+1)-2
       %%
       for in=kn:(((max_value+1)-3)+kn)
           %%
           set_node_lefts(total)=(((max_value+1)^2)+2)+in;
           fprintf(inpFile,'\n*Nset, nset=LeftS-%d, instance=DREAM-1\n',total);
           fprintf(inpFile,'%d,\n',set_node_lefts(total)');
           total=total+1;
       end
       total_all=total_all+total;
       kn=kn+(max_value+1)^2;
   end

   %right face x-axis 
   kn=2;
   total=1;
   total_all=1;
   %back right edge
    for in=kn:(((max_value+1)-3)+kn)
        %%
        set_node_backright(total)=in;
        fprintf(inpFile,'\n*Nset, nset=B_R_edge-%d, instance=DREAM-1\n',total);
        fprintf(inpFile,'%d,\n',set_node_backright(total)');
        total=total+1;
    end
    total=1;
    kn=kn+(max_value+1);
       
    %Back surface
    for jn=1:(max_value-1)
        %%
       for in=kn:(((max_value+1)-3)+kn)
           %%
           back_surface(total)=in;
           fprintf(inpFile,'\n*Nset, nset=BackS-%d, instance=DREAM-1\n',total);
           fprintf(inpFile,'%d,\n',back_surface(total)');
           total=total+1;
       end
              
       total_all=total_all+total;
       kn=kn+(max_value+1);
   end
   total=1;
   
   %back left edge
   for in=kn:(((max_value+1)-3)+kn)
       %%
       set_node_backleft(total)=in;
       fprintf(inpFile,'\n*Nset, nset=B_L_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_backleft(total)');
       total=total+1;
   end
    
   %left face x-axis
   kn=((max_value+1)^2)*(max_value)+2;
   total=1;
   total_all=1;
   
   %Front right edge
   for in=kn:(((max_value+1)-3)+kn)
       %%
       set_node_frontright(total)=in;
       fprintf(inpFile,'\n*Nset, nset=F_R_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_frontright(total)');
       total=total+1;
   end
   total_all=total_all+total;
   kn=kn+(max_value+1);
   
   %front surface
   total=1;
   for jn=1:(max_value-1)
       %%
       for in=kn:(((max_value+1)-3)+kn)
           %%
           set_node_fronts(total)=in;
           fprintf(inpFile,'\n*Nset, nset=FrontS-%d, instance=DREAM-1\n',total);
           fprintf(inpFile,'%d,\n',set_node_fronts(total)');
           total=total+1;
       end
       total_all=total_all+total;
       kn=kn+(max_value+1);
   end
   
   %front left edge
   total=1;
   for in=kn:(((max_value+1)-3)+kn)
       %%
       set_node_frontleft(total)=in;
       fprintf(inpFile,'\n*Nset, nset=F_L_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_frontleft(total)');
       total=total+1;
   end
   total_all=total_all+total;
      
   %corner 8
   kn=1;
   total=1;
   total_all=1;
   for in=kn:kn
       %%
       set_node_C8(total)=in;
       fprintf(inpFile,'\n*Nset, nset=C8-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d\n',set_node_C8(total)');
       total=total+1;
   end
   total_all=total_all+total;
   kn=kn+(max_value+1);
   total=1;
   
   %back bottom edge
   for in=kn:(max_value+1):(max_value)*(max_value+1)
       %%
       set_node_bbottomedge(total)=in;
       fprintf(inpFile,'\n*Nset, nset=B_B_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d\n',set_node_bbottomedge(total)');
       total=total+1;
   end
   total_all=total_all+total;
   kn=kn+(max_value)*(max_value+1)-(max_value+1);
       
   %corner 5
   total=1;
   for in=kn:kn
       %%
       set_node_C5(total)=in;
       fprintf(inpFile,'\n*Nset, nset=C5-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d\n',set_node_C5(total)');
       total=total+1;
   end
   
   total_all=total_all+total;
   kn=kn+(max_value+1);
   total=1;  
   %right bottom edge
   for jn=1:(max_value-1)
       %%
       for in=kn:(max_value)*(max_value+1)+kn:(max_value)*(max_value+1)+kn
           %%
           set_node_rbottomedge(total)=in;
           fprintf(inpFile,'\n*Nset, nset=R_B_edge-%d, instance=DREAM-1\n',total);
           fprintf(inpFile,'%d\n',set_node_rbottomedge(total)');
           total=total+1;
       end
       total_all=total_all+total;
       kn=kn+(max_value+1)^2;
   end
   
   total=1;
   %C7
   for in=kn:kn
       %%
       set_node_C7(total)=in;
       fprintf(inpFile,'\n*Nset, nset=C7-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d\n',set_node_C7(total)');
       total=total+1;
   end
   total_all=total_all+total;
   kn=kn+max_value+1;

   %front bottom edge
   total=1;
   for in=kn:(max_value+1):(max_value)*(max_value+1)+kn-2*(max_value+1)
       %%
       set_node_F_B_edge(total)=in;
       fprintf(inpFile,'\n*Nset, nset=F_B_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d\n',set_node_F_B_edge(total)');
       total=total+1;
   end
   total_all=total_all+total;
   kn=1+((max_value)*(max_value+1))+(max_value+1)^2;
   
   %left bottom edge
   total=1;
   for in=kn:(max_value+1)^2:((max_value+1)^2)*(max_value)
       %%
       set_node_L_B_edge(total)=in;
       fprintf(inpFile,'\n*Nset, nset=L_B_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d\n',set_node_L_B_edge(total)');
       total=total+1;
   end
   total_all=total_all+total;
   kn=kn+(max_value+1)^2;
 
   %bottom surface    
   kn=1;
   kn=kn+(max_value+1)+(max_value+1)^2;
   total=1;
   for jn=1:(max_value-1)
       %%
       for in=kn:(max_value+1)^2:((max_value+1)^2)*(max_value)
           %%
           set_node_bottom(total)=in;
           fprintf(inpFile,'\n*Nset, nset=BotS-%d, instance=DREAM-1\n',total);
           fprintf(inpFile,'%d\n',set_node_bottom(total)');
           total=total+1;
       end
       total_all=total_all+total;
       kn=kn+(max_value+1);
   end
   
   %C6
   total=1;
   kn=kn+((max_value+1)^2)*(max_value-1);
   for in=kn:kn
       %%
       set_node_C6(total)=in;
       fprintf(inpFile,'\n*Nset, nset=C6-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d\n',set_node_C6(total)');
       total=total+1;
   end
   total_all=total_all+total;
       
   %top y-axis 
   %C4
   kn=(max_value+1);
   total=1;
   total_all=1;
   for in=kn:kn
       %%
       set_node_C4(total)=in;
       fprintf(inpFile,'\n*Nset, nset=C4-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_C4(total)');
       total=total+1;
   end
   total_all=total_all+total;
   kn=kn+(max_value+1);
 
   
   %back top edge
   total=1;
   for in=kn:(max_value+1):(max_value)*(max_value+1)
       %%
       set_node_B_T_edge(total)=in;
       fprintf(inpFile,'\n*Nset, nset=B_T_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_B_T_edge(total)');
       total=total+1;
   end
   total_all=total_all+total;
       
   %C1
   total=1;
   kn=(max_value+1)^2;
   for in=kn:kn
       %%
       set_node_C1(total)=in;
       fprintf(inpFile,'\n*Nset, nset=C1-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_C1(total)');
       total=total+1;
   end
   total_all=total_all+total;
             
   %right top edge
   total=1;
   kn=kn+(max_value+1);
   for in=kn:(max_value+1)^2:((max_value)^2)*(max_value+1)
       %%
       set_node_R_T_edge(total)=in;
       fprintf(inpFile,'\n*Nset, nset=R_T_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_R_T_edge(total)');
       total=total+1;
   end
   total_all=total_all+total;
       
   %C3 
   kn=kn+((max_value+1)^2)*(max_value-1);
   total=1;
   for in=kn:kn
       %%
       set_node_C3(total)=in;
       fprintf(inpFile,'\n*Nset, nset=C3-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_C3(total)');
       total=total+1;
   end
   total_all=total_all+total;
       
   %front top edge
   total=1;
   kn=kn+(max_value+1);
   for in=kn:(max_value+1):(max_value)*(max_value-1)+kn
       %%
       set_node_F_T_edge(total)=in;
       fprintf(inpFile,'\n*Nset, nset=F_T_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_F_T_edge(total)');
       total=total+1;
   end
   total_all=total_all+total;
       
   kn=kn+(max_value+1)*(max_value-1);
   total=1;
   %C2
   for in=kn:kn
       %%
       set_node_C2(total)=in;
       fprintf(inpFile,'\n*Nset, nset=C2-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_C2(total)');
       total=total+1;
   end
   total_all=total_all+total;
       
   kn=((max_value+1)^2)+(max_value+1)^2;
   total=1;
   %left top edge
   for in=kn:(max_value+1)^2:((max_value+1)^2)*(max_value-2)+kn
       %%
       set_node_L_T_edge(total)=in;
       fprintf(inpFile,'\n*Nset, nset=L_T_edge-%d, instance=DREAM-1\n',total);
       fprintf(inpFile,'%d,\n',set_node_L_T_edge(total)');
       total=total+1;
   end
   total_all=total_all+total;
      
   %top surface
   kn=(max_value+1)^2 + (max_value+1)*2;
   total=1;  
   for jn=1:(max_value-1)
       %%
       for in=kn:(max_value+1):(max_value)*(max_value-1)+kn
           %%
           set_node_TopS(total)=in;
           fprintf(inpFile,'\n*Nset, nset=TopS-%d, instance=DREAM-1\n',total);
           fprintf(inpFile,'%d,\n',set_node_TopS(total)');
           total=total+1;
       end
       total_all=total_all+total;
       kn=kn+(max_value+1)^2;
    end

    % Applying periodic boundary equations for corners, edges, and
    % surfaces.  The method used to apply the equations is based on the
    % work in:
    
    % Li, S. and A. Wongsto (2004)
    % "Unit cells for micromechanical analyses of particle-reinforced composites." 
    % Mechanics of Materials 36(7),543-572.

    %% UC-UD=FCD
    for node=1:length(set_node_TopS)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nTopS-%d, 1, 1',node);
        fprintf(inpFile,'\nBotS-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodey, 1, 1');
    end
    for node=1:length(set_node_TopS)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nTopS-%d, 2, 1',node);
        fprintf(inpFile,'\nBotS-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
    end
     for node=1:length(set_node_TopS)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nTopS-%d, 3, 1',node);
        fprintf(inpFile,'\nBotS-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodey, 3, 1');
     end
     
    %% UA-UB=FAB 
   
    for node=1:length(set_node_fronts)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nFrontS-%d, 1, 1',node);
        fprintf(inpFile,'\nBackS-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
    end
    
    for node=1:length(set_node_fronts)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nFrontS-%d, 2, 1',node);
        fprintf(inpFile,'\nBackS-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodex, 2, 1');
    end
    
    for node=1:length(set_node_fronts)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nFrontS-%d, 3, 1',node);
        fprintf(inpFile,'\nBackS-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodex, 3, 1');
    end
         
    %% UE-UF=FEF
    for node=1:length(set_node_lefts)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nLeftS-%d, 1, 1',node);
        fprintf(inpFile,'\nRightS-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodez, 1, 1');
    end
    for node=1:length(set_node_lefts)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nLeftS-%d, 2, 1',node);
        fprintf(inpFile,'\nRightS-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodez, 2, 1');
    end
    for node=1:length(set_node_lefts)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nLeftS-%d, 3, 1',node);
        fprintf(inpFile,'\nRightS-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end
      
    %% UIII-UI=FAB+FCD
    for node=1:length(set_node_F_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_T_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
    end
   
    for node=1:length(set_node_F_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_T_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
    end
   
    for node=1:length(set_node_F_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_T_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodex, 3, 1');
    end
   
   %% UII-UI=FAB
   for node=1:length(set_node_B_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_B_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
   end     
   for node=1:length(set_node_B_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_B_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodex, 2, 1');    
   end
   for node=1:length(set_node_B_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_B_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodex, 3, 1');
   end

    %% UIV-UI=FCD

    for node=1:length(set_node_B_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nB_T_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodey, 1, 1');
    end        
    for node=1:length(set_node_B_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nB_T_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
        
    end
    for node=1:length(set_node_B_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nB_T_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nB_B_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodey, 3, 1');    
    end

    %% UVIII-UV=FEF
    for node=1:length(set_node_frontleft)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nB_L_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodez, 1, 1');
    end
    for node=1:length(set_node_frontleft)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nB_L_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodez, 2, 1');
    end
    for node=1:length(set_node_frontleft)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nB_L_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end

    %% UVII-UV=FAB+FEF
    for node=1:length(set_node_frontleft)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_L_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
    end
    for node=1:length(set_node_frontleft)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_L_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodez, 2, 1');
    end
    for node=1:length(set_node_frontleft)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_L_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end

    %% UVI-UV=FAB
    for node=1:length(set_node_backright)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_R_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
    end
    
    for node=1:length(set_node_backright)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_R_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodex, 2, 1');
    end
    
    for node=1:length(set_node_backright)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nF_R_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nB_R_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodex, 3, 1');
    end

    %% UX-UIX=FCD
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nR_T_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodey, 1, 1');
    end
    
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nR_T_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
    end
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nR_T_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodey, 3, 1');
    end
    
    %% UXI-UIX=FCD+FEF
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nL_T_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodey, 1, 1');
    end
   
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nL_T_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
    end
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nL_T_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end

    %% UXII-UIX=FEF
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nL_B_edge-%d, 1, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodez, 1, 1');
    end
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nL_B_edge-%d, 2, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodez, 2, 1');
    end
    for node=1:length(set_node_L_T_edge)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nL_B_edge-%d, 3, 1',node);
        fprintf(inpFile,'\nR_B_edge-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end

    %% U2-U1=FAB     
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC7-%d, 1, 1',node);
        fprintf(inpFile,'\nC8-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC7-%d, 2, 1',node);
        fprintf(inpFile,'\nC8-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodex, 2, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC7-%d, 3, 1',node);
        fprintf(inpFile,'\nC8-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodex, 3, 1');
    end

    %% U3-U1=FAB+FCD
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC3-%d, 1, 1',node);
        fprintf(inpFile,'\nC8-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC3-%d, 2, 1',node);
        fprintf(inpFile,'\nC8-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC3-%d, 3, 1',node);
        fprintf(inpFile,'\nC8-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodex, 3, 1');
    end
  
    %% U4-U1=FCD
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC4-%d, 1, 1',node);
        fprintf(inpFile,'\nC8-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodey, 1, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC4-%d, 2, 1',node);
        fprintf(inpFile,'\nC8-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC4-%d, 3, 1',node);
        fprintf(inpFile,'\nC8-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodey, 3, 1');
    end
    %% U5-U1=FEF
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC5-%d, 1, 1',node);
        fprintf(inpFile,'\nC8-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodez, 1, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC5-%d, 2, 1',node);
        fprintf(inpFile,'\nC8-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodez, 2, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC5-%d, 3, 1',node);
        fprintf(inpFile,'\nC8-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end
    
    %% U6-U1=FAB+FEF
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC6-%d, 1, 1',node);
        fprintf(inpFile,'\nC8-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
    end
  
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC6-%d, 2, 1',node);
        fprintf(inpFile,'\nC8-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodez, 2, 1');
    end
  
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC6-%d, 3, 1',node);
        fprintf(inpFile,'\nC8-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end

    %% U7-U1=FAB+FCD+FEF
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC2-%d, 1, 1',node);
        fprintf(inpFile,'\nC8-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodex, 1, 1');
    end
  
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC2-%d, 2, 1',node);
        fprintf(inpFile,'\nC8-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
    end
  
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC2-%d, 3, 1',node);
        fprintf(inpFile,'\nC8-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end

    %% U8-U1=FCD+FEF
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC1-%d, 1, 1',node);
        fprintf(inpFile,'\nC8-%d, 1, -1',node);
        fprintf(inpFile,'\nrefnodey, 1, 1');
    end
    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC1-%d, 2, 1',node);
        fprintf(inpFile,'\nC8-%d, 2, -1',node);
        fprintf(inpFile,'\nrefnodey, 2, 1');
    end

    for node=1:length(set_node_C5)
        fprintf(inpFile,'\n**Constraint: Constraint-%d\n*Equation',node);
        fprintf(inpFile,'\n3');
        fprintf(inpFile,'\nC1-%d, 3, 1',node);
        fprintf(inpFile,'\nC8-%d, 3, -1',node);
        fprintf(inpFile,'\nrefnodez, 3, 1');
    end
    
    %% Closing the assembly component of the input file
    
    fprintf(inpFile,'\n*End Assembly');
    fprintf(inpFile, '\n**MATERIALS\n**');
    
    %import the centroid for each grain
    centroid=xlsread('inputfile_info.xlsx','centroid');
     
    %import material parameters to be used in the development of materials
    %for each grain.
    [A]=xlsread('inputfile_info.xlsx','Material_parameters');
    
    %% Finalising the input file
%*******************IMPORTANT NOTE*******************************
% The value for 'Depvar' is currently calculated using 12*NSPTL+2,
%where NSPTL is the total number of slip systems across all slip sets
%****************************************************************

   
    %updating the material parameters with the defined local vectors.
    A(57:59)=v1;
    A(65:67)=v2;
    for ii=1:length(grain_order) 
        fprintf(inpFile, '\n*Material, name=MATERIAL-GRAIN%d',ii);
        fprintf(inpFile, '\n*Depvar\n10000,');
        fprintf(inpFile, '\n*User Material, constants=175\n');
        %updating the material parameters with global vectors
        A(60:62)=[u(grain_order(ii)),v(grain_order(ii)),w(grain_order(ii))];
        A(68:70)=[h(grain_order(ii)),k(grain_order(ii)),l(grain_order(ii))];
        %adding euler angles in radians to be used in the UMAT to calculate
        %the angle of the grain with respect to the global system
        A(169:171)=[deg2rad(euler_angle1(ii)),deg2rad(euler_angle2(ii)),deg2rad(euler_angle3(ii))];
        %adding the centroid information in x,y,z coordinates
        A(172:174)=centroid(grain_order(ii),:);
        %adding the calculated equivalent spherical diameter for each grain
        A(175)=diameter(grain_order(ii));
        %printing this information to file
        fprintf(inpFile, '%u, %u, %u, %u, %u, %u, %u, %u\n',A);
    
    end

    fprintf(inpFile,'\n**');
    fprintf(inpFile, '\n**\n** STEP: Loading\n**\n*Step, name=Loading, inc=10000\n*Static\n0.01, 150., 1e-05, 1.');
    fprintf(inpFile, '\n**\n** OUTPUT REQUESTS\n**');
    fprintf(inpFile, '\n*Restart, write, frequency=0\n**');
    fprintf(inpFile, '\n** FIELD OUTPUT: F-Output-1\n**\n*Output, field, variable=PRESELECT\n**');
    fprintf(inpFile, '\n** HISTORY OUTPUT: H-Output-1\n**\n*Output, history, variable=PRESELECT\n**');
    fprintf(inpFile, '\n*End Step');
    
    % close the file
    fclose(inpFile);  
   
    
    toc
   end