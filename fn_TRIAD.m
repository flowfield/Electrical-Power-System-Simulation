function [ result_DCM ] = fn_TRIAD( Reference_Vector, Body_Vector )

V1 = Reference_Vector(:,1) ;
V2 = Reference_Vector(:,2) ;
W1 = Body_Vector(:,1) ;
W2 = Body_Vector(:,2) ;

r1 = V1 ;
temp = cross(V1,V2) ;
r2 = temp / norm(temp) ;
r3 = cross(V1,temp) / norm(temp) ;

s1 = W1 ;
temp = cross(W1,W2) ;
s2 = temp / norm(temp) ;
s3 = cross(W1,temp) / norm(temp) ;

TRIAD_DCM = s1*r1' + s2*r2' + s3*r3' ;

result_DCM = TRIAD_DCM ;
