function [Abs] = FresnelAbs( n2, beta)
% Calculate the reflected fraction of light, according to the Fresnel
% equations
% usage: [R] = FresnelReflec(n1, n2, cosi)
%       n1 : refractive index of 'from' material
%       n2 : refractive index of 'to' material
%     cosi : Cosine of angle of incidence, i.e. dot product between the
%             incident direction and normal vector.
%             This may be a vector or matrix. Values should be positive!
%        R : reflected fraction of light, in the same shape as cosi
%
% Note: Polarization is being ignored here, as well as total internal
% reflection
%


n1=1;

% modified on 16 Jan 2020
beta=90-beta;

% n2=1.54;
%  figure;
% hold on;

% for n2=1.5:0.2:3.2

% for theta=0:90
   cosi=cosd(beta);

cost = sqrt(1 - (n1/n2)^2*(1 - cosi.^2));
Rs = ((n1*cosi - n2*cost) ./ (n1*cosi + n2*cost)).^2;
Rp = ((n1*cost - n2*cosi) ./ (n1*cost + n2*cosi)).^2;
R = (Rs + Rp)/2;

% plot(90-beta,1-R,'k*')

Abs=1-R;
% end

% end