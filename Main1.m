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