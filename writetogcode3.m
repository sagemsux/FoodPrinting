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