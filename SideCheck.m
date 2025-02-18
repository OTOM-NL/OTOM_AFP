function [Z] = SideCheck(p1,p2,p3,P)
Z=0;
cp11 = cross(p2-p1,P-p1);
cp12 = cross(p2-p1,p3-p1);
cp21 = cross(p3-p2,P-p2);
cp22 = cross(p3-p2,p1-p2);
cp31 = cross(p1-p3,P-p3);
cp32 = cross(p1-p3,p2-p3);
d1 = dot(cp11,cp12);
d2 = dot(cp21,cp22);
d3 = dot(cp31,cp32);
if d1>=0 && d2>=0 && d3>=0;
    Z=1;
end;