#include "discequations.cpp"
#include "allocate.cpp"



int main()
{




double residual=1;
int index;

bool loop_switch=true;





int t=0;
while (loop_switch)
	{


        solve_x_momentum_equation(a_W_u, a_E_u, a_N_u , a_S_u , d_p_u, a_p_u, ustar, vstar, pstar);
        solve_y_momentum_equation(a_W_v, a_E_v, a_N_v, a_S_v, d_p_v ,a_p_v, ustar , vstar, pstar);
        solve_pressure_correction_equation(a_E_p , a_W_p , a_N_p , a_S_p , d_p_v , d_p_u , a_p_p , b, ustar,vstar,pprime,niter);
        compute_pressure_velocity_corrections(ustar,u,vstar,v,p,pstar, d_p_u, d_p_v ,pprime);
        compute_boundary_conditions(u,ustar,v,vstar,p,pstar);
        compute_residual(b,residual);
        swap_variables(u,ustar,v,vstar,pstar,p);

		if(residual<tol){
		        std::cout<<"converged to a tolerance of: "<<tol<<" in "	<<t<<" Iterations"<<std::endl;		
			write_file(u,v,p,nx,ny,t);
			 index=t;  
			std::cout<<" "<<std::endl;
			if(t!=0){break;}


			 }

		if(t%interval==0){
		std::cout<<"|| iteration number is: "<<t<<" || continuity residual is "<<residual<< " || "<<std::endl; 
		}
t++;




    }


return 0;

}
