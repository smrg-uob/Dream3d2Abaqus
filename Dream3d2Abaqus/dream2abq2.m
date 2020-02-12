function dream2abq2(voxFileName2, voxFileName, inpFileName)
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
        % - the nodes belonging to each feature and the corresponding
            % feature that node is shared with.  This information is outputed
            % to text files for the nodes x,y, and z component.
    
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

    %open grain boundary distance file
    fidnew=fopen(voxFileName2,'rt');
    rawData2 = textscan(fidnew, '%f,%f');
    fclose(fidnew);

    % load euler angles, coordinates, grain and phase IDs
    euler   = cell2mat(rawData(1:3));
    xyz     = cell2mat(rawData(4:6));
    grains  = cell2mat(rawData(7));
    phases  = cell2mat(rawData(8));
    gbdist  = cell2mat(rawData2(1:2));
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
    %% Opening a file to store element id with gb distance
    gbfile=fopen('gb_element.txt','wt');
    elfile=fopen('el_element.txt','wt');
    
    

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
    %initialising val array
    for i=1:length(grains)
        val(i,1)=1;
    end
    for ii = 1:numel(unique(grains))
        %%
        fprintf(inpFile,'\n*Elset, elset=GRAIN-%d\n',grain_order(ii));
        fprintf(inpFile,'%d, %d, %d, %d, %d, %d, %d, %d, %d\n',elem(grains==grain_order(ii))');
        numels=0;
        
        fprintf(gbfile,'%d,%d\n',[elem(grains==grain_order(ii)),gbdist(grains==grain_order(ii),1)]');
       
        
        for tt=1:length(elem(grains==grain_order(ii)))
            %%
            numels=numels+1;
        end
       numels_total(grain_order(ii))=numels;
    end
   %% list the nodes corresponding to elements for each feature which are at the grain boundary
   for s=1:length(grain_order)
    for i=1:length(grain_order)
        fe(i,1,s)=1;
        %bnd_feature(i,1)=1;
    end
   end
   
    for i=1:length(grain_order)
        elc(i,1)=1;
    end
    
      for ii = 1:numel(unique(grains))
        for t=1:length(elem(:,1))
            if(grains(t)==grain_order(ii))
                %if(gbdist(t,1)==0 && gbdist(t+1,2)~=gbdist(t,2))
                if(gbdist(t,1)==0)
                
                        %bnd_feature(gbdist(t+1,2),fe(gbdist(t+1,2),1,grain_order(ii)),grain_order(ii))=elem(t,1);
                        %fe(gbdist(t+1,2),1,grain_order(ii))=fe(gbdist(t+1,2),1,grain_order(ii))+1;
                        bnd_elements(grain_order(ii),elc(grain_order(ii),1))=elem(t,1);
                        elc(grain_order(ii),1)=elc(grain_order(ii),1)+1;
                    
                  %for i=2:length(elem(1,:))
                    
                   % nodex(grain_order(ii),val(grain_order(ii),1))=coord(elem(t,i),3);
                    %nodey(grain_order(ii),val(grain_order(ii),1))=coord(elem(t,i),2);
                   % nodez(grain_order(ii),val(grain_order(ii),1))=coord(elem(t,i),1);
                    
                    
                    
                    %val(grain_order(ii),1)=val(grain_order(ii),1)+1;
                  %end
                end
             end
            
        end
      end
      
     %% Developing the nodes on the surface of feature with each node defined by which feature it is associated with it
      %organise the nodes on the surface of the feature with corresponding
      %feature these nodes are shared with.
     cc=ones(1000);
      for t=1:length(bnd_elements(:,1))
        for i=1:length(bnd_elements(t,:))
            if((bnd_elements(t,i)) ~= 0)
                for s=1:length(elem(:,1))
                    if(elem(s,1)==bnd_elements(t,i))
                        carraysum=sum(((elem==elem(s,2))|(elem==elem(s,3))|(elem==elem(s,4))|(elem==elem(s,5))|(elem==elem(s,6))|(elem==elem(s,7))|(elem==elem(s,8))|(elem==elem(s,9))),2);
                        carray=((elem==elem(s,2))|(elem==elem(s,3))|(elem==elem(s,4))|(elem==elem(s,5))|(elem==elem(s,6))|(elem==elem(s,7))|(elem==elem(s,8))|(elem==elem(s,9)));
                       
                        %carray=((elem(:,2:9)==elem(s,2))|(elem(:,2:9)==elem(s,3))|(elem(:,2:9)==elem(s,4))|(elem(:,2:9)==elem(s,5))|(elem(:,2:9)==elem(s,6))|(elem(:,2:9)==elem(s,7))|(elem(:,2:9)==elem(s,8))|(elem(:,2:9)==elem(s,9)));
                        totalcarray=[carray,carraysum];
                        fcarrayr=find(totalcarray(:,2:9)>0 & totalcarray(:,10)<8);
                        [fcarrayrr,fcarraycc]=find(totalcarray(:,2:9)>0 & totalcarray(:,10)<8);
                      
                        eles_node=[fcarrayrr,elem(fcarrayr+total_els)];
                        
                        
                        %bndels{t}(feature,cont)=%need to add here the nodes for each feature for this element
                       %checkarray=find(sum(ismember(bnd_elements,fcarrayr),2)>=1);
                       remove_row=bnd_elements;
                       remove_row(t,:)=[0];
                       [checkarrayr,checkarrayc]=find(ismember(remove_row,fcarrayrr')>=1);
                       %checkarray(t)=[];
                       
                       for ls=1:length(checkarrayr)
                           el_val=remove_row(checkarrayr(ls),checkarrayc(ls));
                           for uy=1:length(eles_node(:,1))
                               
                               if(eles_node(uy,1)==el_val)  
                                feature_node(checkarrayr(ls),cc(checkarrayr(ls),t),t)=eles_node(uy,2);
                                cc(checkarrayr(ls),t)=cc(checkarrayr(ls),t)+1;
                               end
                           end
                       end    
                    end
                end
            end        
        end
      end
      
      %% convert 3d array to 2d array for exporting
      %check the size of the array
      feature_node_size=size(feature_node);
      feature_node_2d=reshape(feature_node,[feature_node_size(1),feature_node_size(2)*feature_node_size(3)]);
      
      for i=1:length(feature_node_2d(:,1))
          for iii=1:length(feature_node_2d(1,:))
              if(feature_node_2d(i,iii)~=0)
                nodex(i,iii)=coord(feature_node_2d(i,iii),3);
                nodey(i,iii)=coord(feature_node_2d(i,iii),2);
                nodez(i,iii)=coord(feature_node_2d(i,iii),1);
              else
                nodex(i,iii)=0;
                nodey(i,iii)=0;
                nodez(i,iii)=0;
              end
          end
      end
      
  
   
    fclose(gbfile);
    fclose(elfile);
    
    %% Writing node x cordinate values for each feature
    nodecoordxfile=fopen('node_coordx.txt','wt');
    [rows cols] = size(nodex);
    x = repmat('%f\t',1,(cols-1));
    fprintf(nodecoordxfile,[x,'%f\n'],nodex');
    fclose(nodecoordxfile);
    
    %% Writing node y cordinate values for each feature
    nodecoordyfile=fopen('node_coordy.txt','wt');
    [rows cols] = size(nodey);
    x = repmat('%f\t',1,(cols-1));
    fprintf(nodecoordyfile,[x,'%f\n'],nodey');
    fclose(nodecoordyfile);
    
    %% Writing node x cordinate values for each feature
    nodecoordzfile=fopen('node_coordz.txt','wt');
    [rows cols] = size(nodez);
    x = repmat('%f\t',1,(cols-1));
    fprintf(nodecoordzfile,[x,'%f\n'],nodez');
    fclose(nodecoordzfile);
    
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
        fprintf(inpFile, '\n*Material, name=MATERIAL-GRAIN%d',grain_order(ii));
        fprintf(inpFile, '\n*Depvar\n460,');
        fprintf(inpFile, '\n*User Material, constants=176\n');
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
        A(176)=grain_order(ii);
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