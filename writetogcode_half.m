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

