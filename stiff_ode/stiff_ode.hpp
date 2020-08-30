double *r8vec_linspace_new ( int n, double a, double b );
void stiff_deriv ( double t, double y[], double dydt[] );
void stiff_euler ( double tspan[2], double y0[1], int n, double t[], double y[] );
void stiff_euler_backward ( double tspan[2], double y0[1], int n, double t[], double y[] );
double *stiff_exact ( int n, double t[] );
void stiff_midpoint ( double tspan[2], double y0[1], int n, double t[], double y[] );
void timestamp ( );
