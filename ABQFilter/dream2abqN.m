function dream2abqN(inpFileName)
%% How to use
%Run from the command line using: dream2abqN('name__') where 'name' is
%the name of the input file. Please also include the double underscore at 
%the end of the name you choose.

%Running this function will create two files:

%'name_material.inp'
%'name_sections.inp'

%% Author
%Dylan Agius
%dylan.agius@bristol.ac.uk

%%
    %import euler angles of grains
    eulerinput=regexprep(inpFileName,'__','_euler.csv');
    eulera = xlsread(eulerinput);
    eulera = eulera(2:size(eulera,1),:);


    %create transformation matrix
    zrot=zeros(3,3,size(eulera,1));
    xrot=zeros(3,3,size(eulera,1));
    zrot2=zeros(3,3,size(eulera,1));
    total_rot=zeros(3,3,size(eulera,1));

    for i=1:size(eulera,1)
        zrot(:,:,i)=[cos(eulera(i,1)), sin(eulera(i,1)), 0; -sin(eulera(i,1)), cos(eulera(i,1)),0; 0,0,1];
        xrot(:,:,i)=[1,0,0;0,cos(eulera(i,2)),sin(eulera(i,2));0,-sin(eulera(i,2)),cos(eulera(i,2))];  
        zrot2(:,:,i)=[cos(eulera(i,3)),sin(eulera(i,3)),0;-sin(eulera(i,3)),cos(eulera(i,3)),0;0,0,1];
        total_rot(:,:,i)=transpose(zrot2(:,:,i)*xrot(:,:,i)*zrot(:,:,i));
    end
 
    %vectors in the local coordinate system
    vecs1=[1;0;0];
    vecs2=[0;1;0];

    %rotating local vectors to global system using the transformation matrix
    %developed from the euler angles representing the orientation of each grain
    for i=1:size(eulera,1)
        rotvec1(:,:,i)=(total_rot(:,:,i)*vecs1);
        rotvec2(:,:,i)=(total_rot(:,:,i)*vecs2);
    end

    %creating an input file with the sections required
    %create name
    sectionsinput=regexprep(inpFileName,'__','__sections.inp');

    inpFiles = fopen(sectionsinput,'wt');
    for i=1:size(eulera,1) 
        fprintf(inpFiles,'**Section: Section_Grain_Mat%d\n*Solid Section, elset=poly%d, material=Grain_Mat%d\n,\n',i,i,i);
    end
    fclose(inpFiles);  
  
    %create name
    materialinput=regexprep(inpFileName,'__','__materials.inp');
    %creating an input file for materials
    inpFile = fopen(materialinput,'wt');

    [A]=xlsread('inputfile_info.xlsx','Material_parameters');
   
   
    
    
    %% Finalising the input file
%*******************IMPORTANT NOTE*******************************
% The value for 'Depvar' is currently calculated using 12*NSPTL+2,
%where NSPTL is the total number of slip systems across all slip sets
%****************************************************************

   
    %updating the material parameters with the defined local vectors.
    A(57:59)=vecs1;
    A(65:67)=vecs2;
    for i=1:size(eulera,1)  
        fprintf(inpFile, '\n*Material, name=Grain_Mat%d',i);
        fprintf(inpFile, '\n*Depvar\n460,');
        fprintf(inpFile, '\n*User Material, constants=175\n');
        %updating the material parameters with global vectors
        A(60:62)=rotvec1(:,:,i);
        A(68:70)=rotvec2(:,:,i);
        A(173)=i;
        %%%%%%%%%%%NOTE%%%%%%%%%%%%%%%%
        % These parameters are now obsolete
        %adding euler angles in radians to be used in the UMAT to calculate
        %the angle of the grain with respect to the global system
        %A(169:171)=[deg2rad(euler_angle1(ii)),deg2rad(euler_angle2(ii)),deg2rad(euler_angle3(ii))];
        %adding the centroid information in x,y,z coordinates
        %A(172:174)=centroid(grain_order(ii),:);
        %adding the calculated equivalent spherical diameter for each grain
        %A(175)=diameter(grain_order(ii));
        %printing this information to file
        fprintf(inpFile, '%u, %u, %u, %u, %u, %u, %u, %u\n',A);
    
    end
    
    % close the file
    fclose(inpFile); 
    
    inpfilecreat=regexprep(inpFileName,'__','__.inp');
    
    %elements file, nodes, and elset file
    elemsfile=regexprep(inpFileName,'__','__elems.inp');
    sectsfile=regexprep(inpFileName,'__','__sections.inp');
    nodesfile=regexprep(inpFileName,'__','__nodes.inp');
    elsetfile=regexprep(inpFileName,'__','__elset.inp');
    
    %creating an input file for materials
    inpfile=fopen(inpfilecreat,'wt');

    
    fprintf(inpfile, '*Include, Input = %s\n',materialinput);
    fprintf(inpfile, '*Include, Input = %s\n',elemsfile);
    fprintf(inpfile, '*Include, Input = %s\n',sectsfile);
    fprintf(inpfile, '*Include, Input = %s\n',nodesfile);
    fprintf(inpfile, '*Include, Input = %s', elsetfile);
    
    
    fclose(inpfile);
    
   
    
    
   end