
//domain size

const double lx=1.0;
const double ly=1.0;

//simulation parameters
const int nx=150;
const int ny=150;

const double dx =double( lx /(nx));
const double dy= double( ly /(ny));

const double Re=1000;	//Reynold's number
const double nu=double(1)/Re;  //kinematic viscosity
const  double tol=0.000001;
//solution control parameters

const double urf_pressure=0.3;
const double urf_uvelocity=1.0;
const double urf_vvelocity=1.0;



const double omega_p =1.5;

int niter=60;

const int interval=50;
