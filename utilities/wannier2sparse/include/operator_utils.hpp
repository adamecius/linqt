#ifndef OPERATOR_UTILS 
#define OPERATOR_UTILS

#include <string>
#include <vector>
#include <complex>
using namespace std;

namespace oputil
{

	struct op_matrix
	{
		typedef std::complex<double> value_t;
		
		op_matrix(const int dim )
		{
			data = std::vector< std::vector< value_t > >(dim);
			for(auto& elem : data )
				elem = std::vector< value_t >(dim,0.0); 
		};
		
		value_t& operator()(const int i,const int j)
		{
			return data[i][j];
		}
		
		inline 
		int size()const { return data.size(); }; 

		std::vector< std::vector< value_t > > data;
	};

	const 
	string AVAIL_OPS[] = {	
							"VX"  ,"VY"  ,"VZ","VXSX","VXSY","VXSZ",
							"VYSX","VYSY","VYSZ","VZSX","VZSY","VZSZ",
							"SX", "SY", "SZ", "PUSX", "PUSY", "PUSZ",
							"PDSX", "PDSY", "PDSZ", "TX","TY", "TZ"	
							}; 


	/**
	 * Check if the string labeling an operator belongs to the available operators. 
	 *
	 * @param[out] is_op  true if it belongs false if not
	 * @param[in]  op string containing the operator label.
	 * @param[in]  i0 index which defines the initial position of avail_op list to search for op.
	 * @param[in]  i1 index which defines the final position of avail_op list to search for op.
	 */
	inline
	bool is_operator_in(string op, const int i0,const int i1)
	{
		bool is_op=false;
		for( int i = i0; i < i1; i++)
			is_op+= (op.compare(AVAIL_OPS[i])==0);
		return is_op;
	};

	/**
	 * Check if the string labeling an operator belongs to any of the defined velocities. 
	 *
	 * @param[out] is_op  true if it belongs false if not
	 * @param[in]  op string containing the operator label.
	 */
	inline
	bool is_velocity(string op)
	{
		return is_operator_in(op,0,3);
	};

	/**
	 * Check if the string labeling an operator belongs to any of the defined spin velocities. 
	 *
	 * @param[out] is_op  true if it belongs false if not
	 * @param[in]  op string containing the operator label.
	 */
	inline
	bool is_spinvelocity(string op)
	{
		return is_operator_in(op,3,12);
	};

	/**
	 * Check if the string labeling an operator belongs to any of the defined spin operators. 
	 *
	 * @param[out] is_op  true if it belongs false if not
	 * @param[in]  op string containing the operator label.
	 */
	inline
	bool is_spin(string op)
	{
		return is_operator_in(op,12,15);
	};

        /**
	 * Check if the string labeling an operator belongs to any of the defined spin projection operators. 
	 *
	 * @param[out] is_op  true if it belongs false if not
	 * @param[in]  op string containing the operator label.
	 */
	inline
	bool is_spinprojection(string op)
	{
		return is_operator_in(op,15,21);
	};

	/**
	 * Check if the string labeling an operator belongs to any of the defined torque operators. 
	 *
	 * @param[out] is_op  true if it belongs false if not
	 * @param[in]  op string containing the operator label.
	 */
	inline
	bool is_torque(string op)
	{
		return is_operator_in(op,21,24);
	};

	/**
	 * Extract the direction of the spin fron the operator label. This function does not check if the operator belong to the spins operators. 
	 *
	 * @param[out] direction  The direction of the spin which can be = x,y,z;
	 * @param[in]  op string containing the operator label.
	 */
	inline
	char spin_direction(string op)
	{
		return op.back();
	};

	/**
	 * Extract the direction of the velocity fron the operator label. This function does not check if the operator belong to a velocity operators. 
	 *
	 * @param[out] direction  The direction of the velocity which can be = 0,1,2;
	 * @param[in]  op string containing the operator label.
	 */
	 inline
	 char velocity_direction(string op)
	 {
		 return op[1];
	 };

         /**
	 * Extract the direction of the projection fron the operator label. This function does not check if the operator belong to a projection operator. 
	 *
	 * @param[out] direction  The direction of the projection which can be = +1,-1;
	 * @param[in]  op string containing the operator label.
	 */
	 inline
	 char projection_direction(string op)
	 {
		 return op[1];
	 };

};	


#endif
