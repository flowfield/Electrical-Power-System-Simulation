
function [ lat, lon, alt ] = fn_xyz2lla_ECEF( xyz )

R_E = 6378.145 * 10^3 ;     % m, radius of the Earth

x = xyz(1) ;
y = xyz(2) ;
z = xyz(3) ;

lat = atan2( z, sqrt(x^2+y^2) ) ;
lon = atan2( y, x ) ;
alt = sqrt( x^2 + y^2 + z^2 ) - R_E ;