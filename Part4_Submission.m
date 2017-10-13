% MATLAB code to generate G-code to be input to the food printer to print a
% flower and leaf structure. Note that to run this code, MATLAB R2016b or
% later is required to support local function calling.
% 
% Course: Digital Manufacturing - Dr. Hod Lipson
% Creators: Jacob Joseph, Rayal Raj Prasad, Brian Jin, Jia Huang
% Department of Mechanical Engineering, Columbia University in the City of New York

%% Flower
close all
clear
global gauge
global syringe

%Constants
gauge = 1.29;
syringe = 16;
vol_syringe = syringe^2/4*pi*70;

%Print Settings
resolution = 1000;                  %number of points per revolution traced
walls = 3;                          %number of walls
offset = 1.4;                       %distance between each wall
scale = 20;                         %shape of flower
s = 0.65;                           %size of base cross section
center = [90 60];                   %position of center of flower
layers = 33;                        %number of layers
z_height = linspace(1.05,0.9,layers);   %height per layer
cutoffR = 3;                            %radius of center hole                            

%Dual Extrusion Setting
dualExtrusion = true;
alternate = 3;                      %switch syringe every "alternate" layer


%Function variables
Pflower_final = [];
length = 0;
layerSwitch = 0;
splitRow = 1;

%Booleans
extruder_flag = 0;
first2Syringes = 1;
findRow = 0;

for i=1:layers
    
    %Generate xyz points for the ith layer
    if mod(i,2)
        P = Flower(center,(i-1)*z_height(i),s,walls,-((i-1)/scale)^2+2,offset,resolution,cutoffR);
    else
        P = flip(Flower(center,(i-1)*z_height(i),s,walls,-((i-1)/scale)^2+2,offset,resolution,cutoffR));
    end
    
    %Calculate the distance between each xyz point, store length in P. Also
    %calculates the volume of material that has been extruded
    P(1,4) = 0;
    for j = 2:size(P,1)
        P(j,4) = pdist([P(j,1:3);P(j-1,1:3)],'Euclidean');
        length = length + P(j,4);
    end
    
    ext_volume = length*gauge;
    
    %Check to see if the volume extruded exceeds the volume of 2 syringes.
    %If the condition is met, set the first2Syringes boolean to 0, which
    %will indicate later to split the g-code into 2 files
    if ext_volume > 2*vol_syringe*.8 && first2Syringes
        layerSwitch = i - mod(i,alternate);
        first2Syringes = 0;
    end
    
    %if dual extrusion is selected, then offset x coordinates by 51 mm
    if(extruder_flag==1 && dualExtrusion)
        P = P - repmat([51 -.8 0 0],size(P,1),1);
    end
    
    Pflower_final = [Pflower_final;P];
    
    %insert row of zeros into Pflower_final after every "alternate" layers.
    %A row of zeros tells the writetogcode function when to switch
    %extruders.
    if (mod(i,alternate)==0 && dualExtrusion)
        Pflower_final = [Pflower_final;zeros(1,size(Pflower_final,2))];
        extruder_flag = 1-extruder_flag;
    end
    
    %Finds row in Pflower_final matrix where the g-code should be split
    %into 2 files. Only runs if first2Syringes is 0
    if first2Syringes == 0 && findRow == 0
        splitRow = size(Pflower_final,1);
        findRow = 1;
    end

end

%Plot the Pflower_final points
figure
plot3(Pflower_final(:,1),Pflower_final(:,2),Pflower_final(:,3))
axis equal
%zlim([25 35])

Pflower_final(1,4) = 0;

%Split the g-code if the first2Syringes is false
if first2Syringes
    writetogcode3(Pflower_final);
else
    writetogcode_half(Pflower_final,splitRow,dualExtrusion)
end
%% Leaf
close all
clear all
global gauge
global syringe

%Constants
gauge = 1.29;
syringe = 16;


%Print Settings
resolution = 1000;                      %number of points per layer
walls = 2;                              %number of walls
offset = 2.5;                           %distance between walls

s = 7;                                  %size of leaf
center = [60 60];                       %center position
z_height = 1.1;                         %height per layer
layers = 4;                             %number of layers

Pleaf_final = [];
for i=1:layers
    
    %Generate xyz points for the ith layer
    if mod(i,2)
        P = Leaf1(center,(i-1)*z_height,s,walls,1,offset,resolution);
    else
        P = flip(Leaf1(center,(i-1)*z_height,s,walls,1,offset,resolution));
    end

    Pleaf_final = [Pleaf_final;P];

end

%plot Pleaf_final points
figure
plot3(Pleaf_final(:,1),Pleaf_final(:,2),Pleaf_final(:,3))

Pleaf_final(1,4) = 0;

%Calculate the distance between each point in Pleaf_final
for i = 2:size(Pleaf_final,1)
    Pleaf_final(i,4) = pdist([Pleaf_final(i,1:3);Pleaf_final(i-1,1:3)],'Euclidean');
end

writetogcode_leaf(Pleaf_final)

function [P] = Leaf1(center,z,s,walls,scale,offset,resolution)

%Leaf shape parameters
theta = linspace(-pi,pi,resolution);
leafThickness = 2;
numberOfPoints = 25;


R = zeros(walls,size(theta,2));

%Leaf equation
R(1,:) = scale*s*(leafThickness+.9*cos(8*theta)).*(1+.1*cos(24*theta)).*(.9+.1*cos(numberOfPoints*theta)).*(1+sin(theta));

for i = 2:walls
    R(i,:) = R(i-1,:)+offset;
end

%plot leaf in polar coordinates
polarplot(theta,R)

%convert from polar coordinates to Cartesian coordinates
P_prime = zeros(size(R,1)*size(R,2),2);
row = 1;
for i = 1:size(R,1)
    for j = 1:size(R,2)
        P_prime(row,:) = [R(i,j)*cos(theta(j)) R(i,j)*sin(theta(j))];
        row = row + 1;
    end
end

P_prime = P_prime + repmat(center,size(P_prime,1),1);
P = [P_prime, repmat(z,size(P_prime,1),1)];     %add z coordinate to P

end

function [P] = Flower(center,z,s,walls,scale,offset,resolution,minR)

theta = linspace(0,2*pi,resolution);

%r0-7 equation for each circular shape in flower
r0=1*scale*s;
r1=offset*walls+(2+cos(theta*6+pi))*scale*s;
r2=2*offset*walls+(4+cos(theta*6))*scale*s;
r3=3*offset*walls+(6+cos(theta*6+pi))*scale*s;
r4=4*offset*walls+(8+cos(theta*6))*scale*s;

%vectors for each circular shape, each row represents a wall
r_0 = zeros(walls,size(theta,2));
r_0(1,:) = r0;
r_1 = zeros(size(r_0));
r_1(1,:) = r1;
r_2 = zeros(size(r_0));
r_2(1,:) = r2;
r_3 = zeros(size(r_0));
r_3(1,:) = r3;
r_4 = zeros(size(r_0));
r_4(1,:) = r4;

for i = 2:walls   
    r_0(i,:) = r_0(i-1,:)+offset;
    r_1(i,:) = r_1(i-1,:)+offset;
    r_2(i,:) = r_2(i-1,:)+offset;
    r_3(i,:) = r_3(i-1,:)+offset;
    r_4(i,:) = r_4(i-1,:)+offset;
    
end

%plot the points in polar coordinates
polarplot(theta,r_0,theta,r_1,theta,r_2,theta,r_3,theta,r_4)

R = [r_0;r_1;r_2;r_3];      %consolidate r's into 1 matrix

%Checks to make sure no point is within an area defined by user-specified
%circle. This check ensures that the flower cross section is not too small
%for the printe resolution

deleteRow = 0;
tooSmall = false;
for i = 1:size(R,1)
    
    if min(R(end,:)) < minR
        tooSmall = true;
        break
    end
    
    if min(R(i,:)) > minR
        deleteRow = i;
        break
    end
end

if deleteRow ~= 0
    R = R(deleteRow:end,:);
end

if tooSmall
    R = zeros(1,size(R,1));
end

%Convert polar coordinates to Cartesian coordinates
P_prime = zeros(size(R,1)*size(R,2),2);
row = 1;
for i = 1:size(R,1)
    for j = 1:size(R,2)
        P_prime(row,:) = [R(i,j)*cos(theta(j)) R(i,j)*sin(theta(j))];
        row = row + 1;
    end
end

P_prime = [P_prime, repmat(z,size(P_prime,1),1)];   %add z coordinates to P_prime
P_prime(:,1:3) = (rotz(180)*(P_prime(:,1:3)'))';    %rotate coordinates so to prevent interference with priming line

P = P_prime + repmat([center 0],size(P_prime,1),1);

end

function E = extrusion(L)
global gauge
global syringe
ratio = (gauge^2)/syringe^2;
% e = ratio.*Length;
flag = 0;
E = zeros(1,size(L,1));
for i = 2:size(L,1)
    e = ratio.*L(i);
%     if((V(i,1)-x0)^2 + (V(i,2)-y0)^2) <= (R/4)^2
%         E(i) = E(i-1);
%         flag = 1;
%     else
        if flag == 1
            E(i) = 3*e+E(i-1);
        else
            E(i) = e+E(i-1);
        end
        flag = 0;
        
%    end
end
end

function writetogcode3(V)
global gauge    % gauge and syringe body diameter are defined as global variables
global syringe
initialE = 1;   % initial extrusion value used to initialize priming

fid = fopen('flower.gcode','w');    % creating new gcode file for the flower
fprintf(fid,'M106 S200 \n');    %turn on fan
fprintf(fid,'T0 \n');
fprintf(fid,'G90 \n');  %absolute coordinate system
fprintf(fid,'G21 \n');  %Set to millimeters
fprintf(fid,'G92 X0 Y0 Z0 E0 \n');  %set current position to zero
fprintf(fid,'G0 X0 Y0 Z0 E0 \n');

L = V(:,4); % extracting the Euclidean distances between points

l = 50; % length of the priming line

ratio = (gauge^2)/syringe^2;
E1 = (initialE + ratio*l);  % counter variable used if only a single extruder is used
E = extrusion(L);   % function to compute extruder values between every point

% horizontal priming path of length 'l'
Xstart = V(1,1)-l;
Ystart = V(1,2);
Zstart = V(1,3);
fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',Xstart,Ystart,Zstart,initialE);

Xend = V(1,1);
Yend = V(1,2);
Zend = V(1,3);
fprintf(fid,'G01 X%f Y%f Z%f E%f \n',Xend,Yend,Zend,E1);

% flags
T0 = 1; % flag used during dual extrusion tool change
flag = 0;   % flag to indicate if a row of zeroes is hit
E2 = 1; % counter variable that accumulates the value of the extruder for dual extrusion

% design
for i = 2:size(V,1) % loop runs from the second point
    
    % encountering a row of zeroes i.e. change extruder
    if V(i,1) == 0 && V(i,2) == 0 && V(i,3) == 0 && V(i,4) == 0
        
        % toggling tool flag
        T0 = (1-T0);

        % writing points for the extruder
        if flag == 1
            fprintf(fid,'G0 X0 Y0 Z%f E%f \n',V((i-1),3),E(i)+E2);
        else
            fprintf(fid,'G0 X0 Y0 Z%f E%f \n',V((i-1),3),E(i)+E1);
        end
        
        % depending on the extruder used, the offset is either included or not
        if T0
            fprintf(fid,'T0 \n');
            if flag == 1
                fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',X-l,Y,V((i-1),3),E(i-1)+E2);
            else
                fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',X-l,Y,V((i-1),3),E(i-1)+E1);
            end
            E2 = E2+E1;
            fprintf(fid,'G01 X%f Y%f Z%f E%f \n',X,Y,V((i-1),3),E(i-1)+E2);
        else
            fprintf(fid,'T1 \n');
            if flag == 1
                fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',X-2*l,Y,V((i-1),3),E(i-1)+E2);
                
            else
                fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',X-2*l,Y,V((i-1),3),E(i-1)+E1);
            end
            E2 = E2+E1;
            fprintf(fid,'G01 X%f Y%f Z%f E%f \n',X-51,Y,V((i-1),3),E(i-1)+E2);
        end

        % flag set to 1 to indicate dual extrusion - will need different extrusion values
        flag = 1;
    % logic comes here if it is a non-zero row
    else
        X = V(i,1);
        Y = V(i,2);
        Z = V(i,3);

        % corresponding statements depending on whether dual extrusion is active or not
        if flag == 1
            fprintf(fid,'G01 X%f Y%f Z%f E%f \n',X,Y,Z,E(i)+E2);
        else
            fprintf(fid,'G01 X%f Y%f Z%f E%f \n',X,Y,Z,E(i)+E1);
        end
    end
end

% printing the end of the gcode file
fprintf(fid,'G01 X%f Y%f Z%f E0 \n',Xstart,Ystart,V(end,3));    % go back to home position
fprintf(fid,'M84 \n');  % disable motors
fclose(fid);    % close file
end

function writetogcode_leaf(V)
global gauge    % gauge and syringe body diameter are defined as global variables
global syringe
initialE = 1;   % initial extrusion value used to initialize priming
fid = fopen('leaf.gcode','w');  % creating new gcode file for the leaf                 
fprintf(fid,'M106 S200 \n');    % turn on fan
fprintf(fid,'T0 \n');
fprintf(fid,'G90 \n');  % absolute coordinate system
fprintf(fid,'G21 \n');  % set to millimeters
fprintf(fid,'G92 X0 Y0 Z0 E0 \n');  % set current position to zero
fprintf(fid,'G0 X0 Y0 Z0 E0 \n');

L = V(:,4); % extracting the Euclidean distances between points

l = 50; % length of the priming line

ratio = (gauge^2)/syringe^2;
E1 = (initialE + ratio*l);  % counter variable used if only a single extruder is used
E = extrusion(L);   % function to compute extruder values between every point

% horizontal priming path of length 'l'
Xstart = V(1,1)-l;
Ystart = V(1,2);
Zstart = V(1,3);
fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',Xstart,Ystart,Zstart,initialE);

Xend = V(1,1);
Yend = V(1,2);
Zend = V(1,3);
fprintf(fid,'G01 X%f Y%f Z%f E%f \n',Xend,Yend,Zend,E1);

% flags
T0 = 1; % flag used during dual extrusion tool change
flag = 0;   % flag to indicate if a row of zeroes is hit
E2 = 1; % counter variable that accumulates the value of the extruder for dual extrusion

% design
for i = 2:size(V,1) % loop runs from the second point
    
    % encountering a row of zeroes i.e. change extruder
    if V(i,1) == 0 && V(i,2) == 0 && V(i,3) == 0 && V(i,4) == 0
        
        % toggling tool flag
        T0 = (1-T0);

        % writing points for the extruder
        if flag == 1
            fprintf(fid,'G0 X0 Y0 Z%f E%f \n',V((i-1),3),E(i)+E2);
        else
            fprintf(fid,'G0 X0 Y0 Z%f E%f \n',V((i-1),3),E(i)+E1);
        end
        
        % depending on the extruder used, the offset is either included or not
        if T0
            fprintf(fid,'T0 \n');
            if flag == 1
                fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',X-l,Y,V((i-1),3),E(i-1)+E2);
            else
                fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',X-l,Y,V((i-1),3),E(i-1)+E1);
            end
            E2 = E2+E1;
            fprintf(fid,'G01 X%f Y%f Z%f E%f \n',X,Y,V((i-1),3),E(i-1)+E2);
        else
            fprintf(fid,'T1 \n');
            if flag == 1
                fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',X-2*l,Y,V((i-1),3),E(i-1)+E2);
                
            else
                fprintf(fid,'G01 X%f Y%f Z%f E%f F400 \n',X-2*l,Y,V((i-1),3),E(i-1)+E1);
            end
            E2 = E2+E1;
            fprintf(fid,'G01 X%f Y%f Z%f E%f \n',X-51,Y,V((i-1),3),E(i-1)+E2);
        end

        % flag set to 1 to indicate dual extrusion - will need different extrusion values
        flag = 1;
    % logic comes here if it is a non-zero row
    else
        X = V(i,1);
        Y = V(i,2);
        Z = V(i,3);

        % corresponding statements depending on whether dual extrusion is active or not
        if flag == 1
            fprintf(fid,'G01 X%f Y%f Z%f E%f \n',X,Y,Z,E(i)+E2);
        else
            fprintf(fid,'G01 X%f Y%f Z%f E%f \n',X,Y,Z,E(i)+E1);
        end
    end
end

% printing the end of the gcode file
fprintf(fid,'G01 X%f Y%f Z%f E0 \n',Xstart,Ystart,V(end,3));    % go back to home position
fprintf(fid,'M84 \n');  % disable motors
fclose(fid);    % close file
end

function writetogcode_half(V,split_row,dualExtrusion)
global gauge    % gauge and syringe body diameter are defined as global variables
global syringe
length = size(V,1); % variable to store the number of points in the structure

% find nearest row of zeros (round down)
while(split_row>0 && dualExtrusion)
    if(sum(V(split_row,:))~= 0)
        split_row = split_row - 1;
    else
        break
    end
end
split_row = split_row - 1;  % subtract 1 to get the row before extruder switch

initialE = 1;   % initial extrusion value used to initialize priming
fid1 = fopen('flower_part1.gcode','w'); % flower code - part 1
fid2 = fopen('flower_part2.gcode','w'); % flower code - part 2

L = V(:,4);

% writing file 1
fprintf(fid1,'M106 S200 \n');                    %turn on fan
fprintf(fid1,'T0 \n');              
fprintf(fid1,'G90 \n');                          %absolute coordinate system
fprintf(fid1,'G21 \n');                          %Set to millimeters
fprintf(fid1,'G92 X0 Y0 Z0 E0 \n');              %set current position to zero
fprintf(fid1,'G0 X0 Y0 Z0 E0 \n');

l = 70; % priming line longer to ensure extruders are outside print while switching
ratio = (gauge^2)/syringe^2;
E1 = (initialE + ratio*l); % storing extruder values
E2 = 1;

E = extrusion(L);

% priming path
Xstart = V(1,1)-l;
Ystart = V(1,2);
Zstart = V(1,3);
fprintf(fid1,'G01 X%f Y%f Z%f E%f F400 \n',Xstart,Ystart,Zstart,initialE);


Xend = V(1,1);
Yend = V(1,2);
Zend = V(1,3);
fprintf(fid1,'G01 X%f Y%f Z%f E%f \n',Xend,Yend,Zend,E1);

% flags
T0 = 1;
flag = 0;
Eoffset = 0;

%design
for i = 2:size(V,1)
    
    if i == split_row
       %write head files
        fprintf(fid2,'M106 S200 \n');                    %turn on fan\
        fprintf(fid2,'G90 \n');                          %absolute coordinate system
        fprintf(fid2,'G21 \n');                          %Set to millimeters
        %ensuring correct extrustion value for second file
        fprintf(fid2,'G92 X%f Y%f Z%f E%f \n',Xstart,Y_1,Z_1,E(i-1)+E2);              %set current position to zero
    end
    
    if sum(V(i,:)) == 0
        
        T0 = (1-T0);
        if i<split_row
            if flag == 1 
                fprintf(fid1,'G0 X0 Y0 Z%f E%f \n',V((i-1),3),E(i)+E2);
            else
                fprintf(fid1,'G0 X0 Y0 Z%f E%f \n',V((i-1),3),E(i)+E1);%-Eoffest removed
            end
        else
            if flag == 1
                fprintf(fid2,'G0 X0 Y0 Z%f E%f \n',V((i-1),3),E(i)+E2);
            else
                fprintf(fid2,'G0 X0 Y0 Z%f E%f \n',V((i-1),3),E(i)+E1);%-Eoffest removed
            end
        end
    
        if i<split_row
            if T0
                fprintf(fid1,'T0 \n');
                if flag == 1
                    fprintf(fid1,'G01 X%f Y%f Z%f E%f F400 \n',0,Y,V((i-1),3),E(i-1)+E2); % was X-l
                else
                    fprintf(fid1,'G01 X%f Y%f Z%f E%f F400 \n',0,Y,V((i-1),3),E(i-1)+E1);
                end
                E2 = E2+E1;
                fprintf(fid1,'G01 X%f Y%f Z%f E%f \n',X,Y,V((i-1),3),E(i-1)+E2);
            else
                fprintf(fid1,'T1 \n');
                if flag == 1
                    fprintf(fid1,'G01 X%f Y%f Z%f E%f F400 \n',0,Y,V((i-1),3),E(i-1)+E2); % was X-2*l

                else
                    fprintf(fid1,'G01 X%f Y%f Z%f E%f F400 \n',0,Y,V((i-1),3),E(i-1)+E1);
                end
                E2 = E2+E1;
                fprintf(fid1,'G01 X%f Y%f Z%f E%f \n',X-51,Y,V((i-1),3),E(i-1)+E2);
            end
        else
            if T0
                fprintf(fid2,'T0 \n');
                if flag == 1
                    fprintf(fid2,'G01 X%f Y%f Z%f E%f F400 \n',0,Y,V((i-1),3),E(i-1)+E2);
                else
                    fprintf(fid2,'G01 X%f Y%f Z%f E%f F400 \n',0,Y,V((i-1),3),E(i-1)+E1);
                end
                E2 = E2+E1;
                fprintf(fid2,'G01 X%f Y%f Z%f E%f \n',X,Y,V((i-1),3),E(i-1)+E2);
            else
                fprintf(fid2,'T1 \n');
                if flag == 1
                    fprintf(fid2,'G01 X%f Y%f Z%f E%f F400 \n',0,Y,V((i-1),3),E(i-1)+E2);

                else
                    fprintf(fid2,'G01 X%f Y%f Z%f E%f F400 \n',0,Y,V((i-1),3),E(i-1)+E1);
                end
                E2 = E2+E1;
                fprintf(fid2,'G01 X%f Y%f Z%f E%f \n',X-51,Y,V((i-1),3),E(i-1)+E2);
            end
        end
        flag = 1;
    else
        X = V(i,1);
        Y = V(i,2);
        Z = V(i,3);
        
        % code to store the last value of X,Y and Z for the second code
        if i<split_row
            X_1 = V(i,1);
            Y_1 = V(i,2);
            Z_1 = V(i,3);
        end
        
        if i<split_row

            if flag == 1
                fprintf(fid1,'G01 X%f Y%f Z%f E%f \n',X,Y,Z,E(i)+E2);
            else
                fprintf(fid1,'G01 X%f Y%f Z%f E%f \n',X,Y,Z,E(i)+E1);
            end
        else
            if flag == 1
                fprintf(fid2,'G01 X%f Y%f Z%f E%f \n',X,Y,Z,E(i)+E2);
            else
                fprintf(fid2,'G01 X%f Y%f Z%f E%f \n',X,Y,Z,E(i)+E1);
            end
        end    
    end
    
    
end

% end of file code
fprintf(fid1,'G01 X%f Y%f Z%f E0 \n',Xstart,Y_1,Z_1);
fprintf(fid1,'M84 \n');
fclose(fid1);
fprintf(fid2,'G01 X%f Y%f Z%f E0 \n',Xstart,Ystart,V(end-1,3));
fprintf(fid2,'M84 \n');
fclose(fid2);

end