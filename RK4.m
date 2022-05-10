function result = RK4( Fcn, X, dt )

k1 = Fcn( X ) * dt ;
k2 = Fcn( X + 0.5*k1 ) * dt ;
k3 = Fcn( X + 0.5*k2 ) * dt ;
k4 = Fcn( X + k3 ) * dt ;

k = ( k1 + 2*k2 + 2*k3 + k4 ) / 6 ;

result = X + k ;