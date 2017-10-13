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