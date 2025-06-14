//#include "matrix.h"
#include "params.cpp"
#include <Eigen/Dense>
#include<iostream>
#include<fstream>

using dmatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;   //the loops are accessed in row major so this is a dirty way to improve the performance without upending the code

void solve_x_momentum_equation(dmatrix &a_E_u,dmatrix &a_W_u,dmatrix &a_N_u,dmatrix &a_S_u,dmatrix &d_p_u,dmatrix &a_p_u, dmatrix &ustar, dmatrix &vstar ,dmatrix &pstar);
void solve_y_momentum_equation(dmatrix &a_E_v,dmatrix &a_W_v,dmatrix &a_N_v,dmatrix &a_S_v,dmatrix &d_p_v,dmatrix &a_p_v, dmatrix &ustar, dmatrix &vstar ,dmatrix &pstar);
void compute_pressure_velocity_corrections(dmatrix &ustar,dmatrix &u, dmatrix &vstar ,dmatrix &v,dmatrix &p,dmatrix &pstar,dmatrix &de,dmatrix &dn, dmatrix &pprime ) ;
void write_file(dmatrix &u,dmatrix &v ,dmatrix &p ,int nx,int ny,int number);
void swap_variables(dmatrix &u ,dmatrix &ustar ,dmatrix &v,dmatrix &vstar ,dmatrix &pstar ,dmatrix &p);






void solve_x_momentum_equation(dmatrix &a_E_u,dmatrix &a_W_u,dmatrix &a_N_u,dmatrix &a_S_u,dmatrix &d_p_u,dmatrix &a_p_u, dmatrix &ustar, dmatrix &vstar ,dmatrix &pstar)
{

double F_e_u;
double F_w_u;
double F_n_u;
double F_s_u;

double D_e_u;
double D_w_u;
double D_n_u;
double D_s_u;
double Ae;



for(int i=1;i<=nx-2;i++)
{
    for(int j=1;j<=(ny-1);j++)
		{

		F_e_u=0.5*dy*(ustar(i+1,j)+ustar(i,j));
		F_w_u=0.5*dy*(ustar(i,j)+ustar(i-1,j));
		F_n_u=0.5*dx*(vstar(i,j)+vstar(i+1,j));
		F_s_u=0.5*dx*(vstar(i+1,j-1)+vstar(i,j-1));


		D_e_u=nu*dy/dx;
		D_w_u=nu*dy/dx;
		D_n_u=nu*dx/dy;
		D_s_u=nu*dx/dy;


		//need to find coefficents as per upwind scheme
		a_E_u(i,j)=std::max(0.0,-F_e_u)+   D_e_u*std::max(0.0,std::pow((1-  (0.1*std::abs(F_e_u))/(D_e_u)  )    ,5)       );
		a_W_u(i,j)=std::max(F_w_u,0.0) +   D_w_u*std::max(0.0,std::pow((1-  (0.1*std::abs(F_w_u))/(D_w_u)  )    ,5)       ); 
		a_N_u(i,j)=std::max(0.0,-F_n_u)+   D_n_u*std::max(0.0,std::pow((1-  (0.1*std::abs(F_n_u))/(D_n_u)  )    ,5)       );
		a_S_u(i,j)=std::max(F_s_u,0.0) +   D_s_u*std::max(0.0,std::pow((1-  (0.1*std::abs(F_s_u))/(D_s_u)  )    ,5)       );

		a_p_u(i,j)=a_W_u(i,j)+a_E_u(i,j)+a_N_u(i,j)+a_S_u(i,j)+(F_e_u-F_w_u)+(F_n_u-F_s_u);

		Ae =dy;

		d_p_u(i,j)=Ae/a_p_u(i,j);

	}

}


	for(int i=1;i<=(nx-2);i++)
	{

		for(int j=1;j<=(ny-1);j++)
		{

		ustar(i,j) =(a_E_u(i,j)*ustar(i+1,j) + a_W_u(i,j)*ustar(i-1,j) + a_N_u(i,j)*ustar(i,j+1) + a_S_u(i,j)*ustar(i,j-1))/a_p_u(i,j) - d_p_u(i,j)*(pstar(i+1,j) - pstar(i,j));



		}

	}


}


void solve_y_momentum_equation(dmatrix &a_E_v,dmatrix &a_W_v,dmatrix &a_N_v,dmatrix &a_S_v,dmatrix &d_p_v,dmatrix &a_p_v, dmatrix &ustar, dmatrix &vstar ,dmatrix &pstar)
{


double F_e_v;
double F_w_v;
double F_n_v;
double F_s_v;

double D_e_v;
double D_w_v;
double D_n_v;
double D_s_v;

double An;



for(int i=1;i<=(nx-1);i++)
{	
	for(int j=1;j<=(ny-2);j++)
	{

		F_e_v=0.5*dy*(ustar(i,j)+ustar(i,j+1));
		F_w_v=0.5*dy*(ustar(i-1,j)+ustar(i-1,j+1));
		F_n_v=0.5*dx*(vstar(i,j)+vstar(i,j+1));
		F_s_v=0.5*dx*(vstar(i,j-1)+vstar(i,j));

		//need to find coefficents as per upwind scheme

		D_e_v=nu*dy/dx;
		D_w_v=nu*dy/dx;
		D_n_v=nu*dx/dy;
		D_s_v=nu*dx/dy;


		a_E_v(i,j)=std::max(0.0,-F_e_v)+   D_e_v*std::max(0.0,std::pow((1-  (0.1*std::abs(F_e_v))/(D_e_v)        )       ,5)            );
		a_W_v(i,j)=std::max(F_w_v,0.0)+    D_w_v*std::max(0.0,std::pow((1-  (0.1*std::abs(F_w_v))/(D_w_v)        )       ,5)            );
		a_N_v(i,j)=std::max(0.0,-F_n_v)+   D_n_v*std::max(0.0,std::pow((1-  (0.1*std::abs(F_n_v))/(D_n_v)        )       ,5)            );
		a_S_v(i,j)=std::max(F_s_v,0.0) +   D_s_v*std::max(0.0,std::pow((1-  (0.1*std::abs(F_s_v))/(D_s_v)        )       ,5)            );

		a_p_v(i,j)=a_W_v(i,j)+a_E_v(i,j)+a_N_v(i,j)+a_S_v(i,j)+(F_e_v-F_w_v)+(F_n_v-F_s_v);

		An=dx;

		d_p_v(i,j)=An/a_p_v(i,j);

	}

}



for(int i=1;i<=(nx-1);i++)
{
	for(int j=1;j<=(ny-2);j++)
	{

	vstar(i,j)=(a_E_v(i,j)*vstar(i+1,j) + a_W_v(i,j)*vstar(i-1,j) + a_N_v(i,j)*vstar(i,j+1) + a_S_v(i,j)*vstar(i,j-1))/a_p_v(i,j) - d_p_v(i,j)*(pstar(i,j+1) - pstar(i,j));



	}
}



}

void solve_pressure_correction_equation(dmatrix &a_E_p,dmatrix &a_W_p,dmatrix &a_N_p,dmatrix &a_S_p,dmatrix &d_p_v, dmatrix &d_p_u, dmatrix &a_p_p, dmatrix &b, dmatrix &ustar, dmatrix &vstar,dmatrix &pprime,int niter)
{


for(int i=1;i<=(nx-1);i++)
{
    	for(int j=1;j<=(ny-1);j++)

		{

		a_E_p(i,j) = d_p_u(i,j)*dy;
        	a_W_p(i,j) = d_p_u(i-1,j)*dy;
        	a_N_p(i,j) = d_p_v(i,j)*dx;
        	a_S_p(i,j)= d_p_v(i,j-1)*dx;

        	a_p_p(i,j) = a_E_p(i,j) + a_W_p(i,j) + a_N_p(i,j) + a_S_p(i,j);
			
		 b(i,j) =-ustar(i,j)*dy + ustar(i-1,j)*dy -vstar(i,j)*dx +vstar(i,j-1)*dx;

	

		}
		
}




for(int i=0;i<nx+1;i++)
{
	for(int j=0;j<ny+1;j++)
	{

		pprime(i,j)=0;

	}

}

//solving for pprime by SOR

for(int iter_count=1;iter_count<niter;iter_count++)
{

#pragma omp parallel for schedule(static)
for(int i=1;i<=(nx-1);i++)
{

	for(int j=1;j<=(ny-1);j++)
	{

		pprime(i,j)=(1-omega_p)*pprime(i,j)+omega_p* (a_E_p(i,j)*pprime(i+1,j) + a_W_p(i,j)*pprime(i-1,j) + a_N_p(i,j)*pprime(i,j+1) + a_S_p(i,j)*pprime(i,j-1) + b(i,j))/a_p_p(i,j);

	}

}

}




}




void compute_pressure_velocity_corrections(dmatrix &ustar,dmatrix &u, dmatrix &vstar ,dmatrix &v,dmatrix &p,dmatrix &pstar,dmatrix &de,dmatrix &dn, dmatrix &pprime ) 
{

for(int i=1;i<=(nx-1);i++)
{
	for(int j=1;j<=(ny-1);j++)
	{
		p(i,j)=pstar(i,j)+urf_pressure*pprime(i,j);

	}
}


for(int i=1;i<=nx-2;i++)
{
    for(int j=1;j<=(ny-1);j++)
{


		u(i,j)=ustar(i,j)+(urf_uvelocity)*de(i,j)*(pprime(i,j)-pprime(i+1,j));



	}

}


for(int i=1;i<=(nx-1);i++)
{	
	for(int j=1;j<=(ny-2);j++)
	{

		v(i,j)= vstar(i,j)+(urf_vvelocity)*dn(i,j)*(pprime(i,j)-pprime(i,j+1));


	}

}


}


void compute_residual(dmatrix &b,double &residual)
{

 residual=0;

for(int i=0;i<=nx-1;i++)
{
	for(int j=0;j<=ny-1;j++)
	{
	
	residual=residual+std::abs(b(i,j));	
	}
}


}




void compute_boundary_conditions(dmatrix &u ,dmatrix &ustar ,dmatrix &v,dmatrix &vstar,dmatrix &p ,dmatrix&pstar)
{

for (int i=0; i<(nx); i++) 
{

    	ustar(i,0) = -ustar(i,1);
    	ustar(i,ny) = 2 - ustar(i,ny-1);


	u(i,0) = -u(i,1);
    	u(i,ny) = 2 - u(i,ny-1);


}


for(int j=0;j<ny+1;j++)
{

	ustar(0,j)=0;
	ustar(nx-1,j)=0;

	u(0,j)=0;
	u(nx-1,j)=0;

}





for (int j=0;j<(ny); j++)
{
    vstar(0,j) =- vstar(1,j);
    vstar(nx,j) =- vstar(nx-1,j);	

    v(0,j) =- v(1,j);
    v(nx,j) =- v(nx-1,j);	


}

for (int i=0; i <(nx+1); i++)
{

	vstar(i,0)=0.0;
	vstar(i,ny-1)=0;

	v(i,0)=0.0;
	v(i,ny-1)=0;	
 		
}



for(int i=0;i<nx+1;i++)
{	
	pstar(i,0)=pstar(i,1);
	pstar(i,ny)=pstar(i,ny-1);
	
	p(i,0)=p(i,1);
	p(i,ny)=p(i,ny-1);
	

}



for(int j=0;j<ny+1;j++)
{
	pstar(0,j)=pstar(1,j);
	pstar(nx,j)=pstar(nx-1,j);
	
	p(0,j)=p(1,j);
	p(nx,j)=p(nx-1,j);
	
	
}


}

void write_file(dmatrix &u,dmatrix &v ,dmatrix &p ,int nx,int ny,int number)
{


    std::string str="final";
	std::string ext=".csv";



	//str=str.append(std::to_string(number));
	str=str.append(ext);

	std::fstream data;

	data.open(str,std::ios::out);


dmatrix ucenter(nx,ny);
dmatrix vcenter(nx,ny);
dmatrix pcenter(nx,ny);

	data<<" x , y , z , u , v , p "<<std::endl;

for(int i=0;i<=nx-1;i++)
{
	for(int j=0;j<=ny-1;j++)
	{

		ucenter(i,j)=0.5*(u(i,j)+u(i,j+1));
		vcenter(i,j)=0.5*(v(i,j)+v(i+1,j));
		pcenter(i,j) = 0.25*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1));

	}
}

			


for (int i=0; i<(nx); i++)
{
	for (int j=0; j<(ny); j++)
	
	{
	    
		
        data<<i*dx<<" , "<<j*dy<<" , "<<"0"<<" , "<<ucenter(i,j)<<" , "<<vcenter(i,j)<<" , "<<pcenter(i,j)<<std::endl;
					
	}
	
}

    data.close();

}


void swap_variables(dmatrix &u ,dmatrix &ustar ,dmatrix &v,dmatrix &vstar ,dmatrix &pstar ,dmatrix &p)
{



for(int i=0;i<(nx);i++)
{

    for(int j=0;j<(ny+1);j++)
	{

		ustar(i,j)=u(i,j);



	}

}



for(int i=0;i<(nx+1);i++)
{

       	for(int j=0;j<(ny);j++)
		{

			vstar(i,j)=v(i,j);



	}

}


for(int i=0;i<(nx+1);i++)
{
	for(int j=0;j<(ny+1);j++)
	{
		
			pstar(i,j)=p(i,j);
			


	}
}



}
