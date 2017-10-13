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