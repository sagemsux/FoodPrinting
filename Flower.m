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