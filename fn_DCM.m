function [ DCM ] = fn_DCM( axis_no, angle )

DCM = eye(3) ;

if ( axis_no == 1 )
    DCM(2,2) = cos(angle) ;
    DCM(2,3) = sin(angle) ;
    DCM(3,2) = - sin(angle) ;
    DCM(3,3) = cos(angle) ;
elseif ( axis_no == 2 )
    DCM(1,1) = cos(angle) ;
    DCM(1,3) = - sin(angle) ;
    DCM(3,1) = sin(angle) ;
    DCM(3,3) = cos(angle) ;
elseif ( axis_no == 3 )
    DCM(1,1) = cos(angle) ;
    DCM(1,2) = sin(angle) ;
    DCM(2,1) = - sin(angle) ;
    DCM(2,2) = cos(angle) ;
end