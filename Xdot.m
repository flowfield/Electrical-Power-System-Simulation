
function [ result ] = Xdot( X )

mu = 398601.2 * (10^3)^3 ;  % m^3/sec^2

Xdot = X(4) ;
Ydot = X(5) ;
Zdot = X(6) ;

r = sqrt( X(1)^2 + X(2)^2 + X(3)^2 ) ;

VXdot = - mu/r^3 * X(1) ;
VYdot = - mu/r^3 * X(2) ;
VZdot = - mu/r^3 * X(3) ;

result = [ Xdot; Ydot; Zdot; VXdot; VYdot; VZdot ] ;