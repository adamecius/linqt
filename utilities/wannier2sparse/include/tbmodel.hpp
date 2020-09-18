#ifndef TBMODELS 
#define TBMODELS
#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <tuple>
#include <complex>
#include <map>
#include <iostream>
#include <limits>
#include <cassert>
#include <functional>
#include "wannier_parser.hpp" //for read_disorder_file
#include "hopping_list.hpp"
#include "operator_utils.hpp"

using namespace std;


/**
 * This functiosn assert that two inputs are equal. the inputs should have overloeaded the operator == . 
 *
 * @param[out] equal.  true if they are both equal, false if not
 * @param[in]  x Input parameter to be compare. It should have the ovearled version of ==.
 * @param[in]  y Input parameter to be compare with X.
 */
template<typename T>
bool assert_equal(const T x,const T y )
{
	bool equal = ( x == y );
    if( !equal )
        std::cerr<<"ASSERT_EQUAL FAILED: "<<x<<" != "<<y<<std::endl;
    assert(x==y);
    return equal;
}


/**
 * Transform a site tag into real space cartesian vector.
 *
 * @param[out] cart_vect. A three-dimensional real vector, constructed from the tag and the lattice vectors
 * @param[in]  tag. An identification of the position of a unit cell. In this version of the code, this is the position in units of lattice vector.
 * @param[in]  lattice_vector. The set of lattice vectors
 */
array<double,3> 
tag2cartesian(const array<int,3>& tag,const array < array<double,3> , 3 >& lat_vecs  );


/**
 * Compute the volume of a given unit cell.
 *
 * @param[out] Volume. the volume of the unit cell
 * @param[in]  UC. Unit cell in the form of an array of arrays.
 */
double 
volume( const array < array<double,3> , 3 >&  uc );


/*! The tbmodel class 
 *  This class is responsible for reading wannier files
 *  and converting it into a hopping_list structure
 *  It is also responsible for constructing other tight-binding operators
 * */

class tbmodel
{
    public:
    typedef tuple<string, array<double, 3> > orbPos_t;
    typedef array < array<double,3> , 3 > unitCell_t;
	typedef multimap< string  ,  hopping_kind > hopping_kind_list ;

    void readUnitCell(const string inputfile);


	void readOrbitalPositions(const string inputfile);

	void readWannierModel(const string inputfile);
    
	void readStaticDisorder(const string inputfile);

	/**
	 * Transform a site tag into real space cartesian vector.
	 *
	 * @param[out] cart_vect. A three-dimensional real vector, constructed from the tag and the lattice vectors
	 * @param[in]  tag. An identification of the position of a unit cell. In this version of the code, this is the position in units of lattice vector.
	 * @param[in]  lattice_vector. The set of lattice vectors
	 */
    map<int,string> map_id2orbs();

	/**
	 * Transform a site tag into real space cartesian vector.
	 *
	 * @param[out] cart_vect. A three-dimensional real vector, constructed from the tag and the lattice vectors
	 * @param[in]  tag. An identification of the position of a unit cell. In this version of the code, this is the position in units of lattice vector.
	 * @param[in]  lattice_vector. The set of lattice vectors
	 */
    map<int,int> map_id2spin();
    
    
    //MAIN FUNCTIONS
	inline
	hopping_list& Hopping_List()
	{
		return hl;
	}


	/**
	 * Extract the local hoppings for the hopping list. 
	 * In this context, the local hoppings are the one satisfying \sum_{i in cell} H_ij|i><j| with 
	 * |i> the basis vector in which the hoppings <i,0,0,0|H|j,0,0,0> are defined. 
	 * The (0,0,0) indexes is to note that this operator only consider hoppings within the unit cell
	 * @return density_hopping_list Returns the density operator in the format of a hopping list.
	 */
    hopping_list createLocalHopping_list( );


	/**
	 * Construct a density operator for the hopping list an a matrix defining an operator. 
	 * In this context, the matrix is defined only in the unit cell, hereby the reason for calling it a density operator.
	 * The density operator is then defined as \sum_{i in cell} M_ij|i><j| with 
	 * |i> the basis vector in which the hoppings <i,0,0,0|H|j,0,0,0> are defined. 
	 * The (0,0,0) indexes is to note that this operator only consider hoppings within the unit cell
	 * @return density_hopping_list Returns the density operator in the format of a hopping list.
	 */
    hopping_list createHoppingDensity_list( oputil::op_matrix OP );

	/**
	 * Construct the current operator for the hopping list. 
	 * In this context, a current operator is defined as \sum_{i,j} |i><j| I*H_{i,j}*(R_i-R_j)_dir with 
	 * |i> the basis vector in which the hoppings <i|H|j> are defined. 
	 * Here i represents both orbital and unit cell indexes.
	 * and d, representes the direction in which you want the current
	 *
	 * @return current_hopping_list Returns the current operator in the format of a hopping list.
	 * @param[in]  Direction. The direction in which you want the current.
	 */
    hopping_list createHoppingCurrents_list(const int dir);

	/**
	 * Construct the density operator in the direction u = ( sin(theta)*cos(phi),sin(theta)*sin(phi), cos(theta); 
	 * In this context,  OP  =   Su, the spin direction.
	 *
	 * @return density_hopping_list Returns the density operator in the format of a hopping list.
	 * @param[in]  Theta. angle defining the altitute direction in which you want the current.
	 * @param[in]  Phi. angle defining the azimuth
	 */
    hopping_list createHoppingSpinDensity_list(const double theta, const double phi);

        /**
	 * Construct the density operator in the direction u = ( sin(theta)*cos(phi),sin(theta)*sin(phi), cos(theta); 
	 * In this context,  OP  =   PsSu, the spin projector direction.
	 *
	 * @return density_hopping_list Returns the density operator in the format of a hopping list.
	 * @param[in]  Theta. angle defining the altitute direction in which you want the current.
	 * @param[in]  Phi. angle defining the azimuth
         * @param[in]  S. eigenvalue of the spin projection.
	 */
  hopping_list createHoppingSpinProjectionDensity_list(const double theta, const double phi, const int s);
    
	/**
	 * Construct the density operator in the direction u = ( sin(theta)*cos(phi),sin(theta)*sin(phi), cos(theta); 
	 * In this context,  OP  =   Su, the spin direction.
	 *
	 * @return density_hopping_list Returns the density operator in the format of a hopping list.
	 * @param[in]  Theta. angle defining the altitute direction in which you want the current.
	 * @param[in]  Phi. angle defining the azimuth
	 */
    hopping_list createHoppingTorqueDensity_list(const double theta, const double phi);

	/**
	 * Construct the current density operator in the direction u = ( sin(theta)*cos(phi),sin(theta)*sin(phi), cos(theta) ); 
	 *
	 * @return spin_current density_hopping_list Returns the Current Density operator in the format of a hopping list.
	 * @param[in]  Theta. angle defining the altitute direction in which you want the current.
	 * @param[in]  Phi. angle defining the azimuth
	 */
    hopping_list createHoppingSpinCurrents_list(const int dir, const double theta, const double phi );


	//DERIVED FUNCTIONS
    hopping_list createHoppingCurrents_list(const std::string op );
    
    hopping_list createHoppingSpinorialDensity_list(const std::array< std::complex<double>,4 > op);

    hopping_list createHoppingSpinDensity_list(const std::string op );

    hopping_list createHoppingSpinProjectionDensity_list(string op);  
   
    hopping_list createHoppingSpinProjectionDensity_list(const int s, const char sdir);

    hopping_list createHoppingTorqueDensity_list(const std::string op );

    hopping_list createHoppingSpinCurrents_list(string op);  
   
    hopping_list createHoppingSpinCurrents_list(const int dir, const char sdir);

    hopping_list WannierOperator(std::string op_id );

	//hopping_list add_onsite_disorder(const hopping_list  hl , double W=0.0 );

	inline 
	void getSpinIndex(const int& n, int& s)const {  s = n/this->NumberOfSites() ; };

	inline 
	int  getSpinIndex(const int& n) const {  return  n/this->NumberOfSites() ; };


	private:
	int NumberOfSites() const { return orbPos_list.size()/MAX_SPIN; } 


	inline
	void ChangeSpinIndex(int &n, const int s) const {  n = n%NumberOfSites() + s*NumberOfSites(); };

	inline
	hopping_list::edge_t ChangeSpinIndex(const hopping_list::edge_t& n,const hopping_list::edge_t& s  ) const 
	{  
		return {  n[0]%NumberOfSites() + s[0]*NumberOfSites(), n[1]%NumberOfSites() + s[1]*NumberOfSites()};
	}
	//inline orbPos_list.size() gives you the number of positions. Including the repeated one due to MAX_SPIN

	oputil::op_matrix createSpinorialOp(const std::array< std::complex<double>,4 > op);

	oputil::op_matrix createSpinMatrix(const double theta, const double phi);

        oputil::op_matrix createSpinProjectionMatrix(const double theta, const double phi, const int s);

	oputil::op_matrix createTorqueMatrix(const double theta, const double phi);

	public:
    int num_orbs;
    unitCell_t lat_vecs;
    vector< orbPos_t > orbPos_list;
    hopping_list hl;

	private :
	const int MAX_SPIN = 2; 
	const int SPIN_UP = 0; 
	const int SPIN_DW = 1; 
	

};

#endif
