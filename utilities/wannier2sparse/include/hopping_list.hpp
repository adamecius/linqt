#ifndef HOPPING_LIST 
#define HOPPING_LIST
#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <tuple>
#include <complex>
#include <map>
#include <fstream>
#include <iostream>
#include <limits>
#include <cassert>
#include <functional>
#include <algorithm>
#include <iostream>
#include <chrono> 
#include "sparse_matrix.hpp"
#include "hopping_kind.hpp"
#include "hopping_list_util.hpp"

using namespace std;
using namespace std::chrono; 

  
/**
 * The hopping_list class handles the hoppings of wannier2sparse.
 * This class parse the hoppings from the input files 
 *
 */
class hopping_list
{
	public: 
	typedef hopping_kind hopkind_t;
    typedef hopkind_t::value_t value_t;
    typedef array<int, 2> edge_t;
    typedef array<int, 3> cellID_t;
    typedef tuple< cellID_t,hopkind_t,edge_t > hopping_t;


	/**
	 * Construct hopping list which posses no wannier functions
	 * and it is assume to exist on a single cell 
	 */
    hopping_list():cellSizes({1,1,1}), num_wann(0){};

	/**
	* Return the Size of the Wannier Basis
	* number.
	*
	* @return real component
	*/
    inline int WannierBasisSize() const
    {
        return this-> num_wann ; 
    };

	/**
	* Return Set the size of the Wannier Basis
	* number.
	*
	* @param[in] num_wann Wannier Basis Size.
	*/
    inline void SetWannierBasisSize(const int num_wann)
    {
        assert( num_wann > 0 );
        this-> num_wann  = num_wann;
        return ; 
    };

	/**
	* Return the number of unit cells needed to define the hoppings. 
	*
	* @return cellSizes A three dimensional array containing the number of unit cells in each spatial direction
	*/   
	cellID_t Bounds()
    {
        return array<int, 3>({cellSizes[0],cellSizes[1],cellSizes[2]});
    };

	/**
	* Set the number of unit cells needed to define the hoppings. 
	*
	* @param[in] cellSizes A three dimensional array containing the number of unit cells in each spatial direction
	*/   
    inline void SetBounds(const cellID_t& cellSizes)
    {
        assert( this-> num_wann > 0 && cellSizes[0]>0&& cellSizes[1]>0&&cellSizes[2]>0 );
        this->cellSizes = cellSizes;
        for(const auto& x: cellSizes ) 
            this->num_wann *=x;
        return ; 
    };


	/**
	* Return the number of unit cells needed to define the hoppings in a given direction. 
	* @param[in] dir The direction index =0,1,2, where the number of unit cells is required.
	* @return cellSize[dir] Number of unit cells in the dir direction.
	*/  
	int cellID_index(const cellID_t cidx )
    {
        const cellID_t bounds = this->Bounds();
       
        int index = 0;
        for( int i = 0; i+1 < cidx.size(); i++ )
            index += ( (cidx[i]+bounds[i])%bounds[i] )*bounds[i+1];
        index += ( cidx.back()+bounds.back() )%( bounds.back() );
        return index;
    };


	/**
	* Add this hopping to the hopping list. 
	* If the hopping already exists, add it as a different hopping with the same tag.
	*
	* @param[in] cellID Destintation cell ID.
	* @param[in] value  Hopping value.
	* @param[in] edge   Initial and final orbital ID.
	*/
	void AddHopping(const cellID_t& _cellID,const value_t& _value ,const edge_t& _edge);

	/**
	* Add this hopping to the hopping list. 
	* If the hopping already exists, add it as a different hopping with the same tag.
	*
	* @param[in] cellID Destintation cell ID.
	* @param[in] value  Hopping Kind.
	* @param[in] edge   Initial and final orbital ID.
	*/
	void AddHopping(const cellID_t& _cellID,const hopping_kind& _hop , const edge_t& _edge);


	/**
	* Add this hopping to the hopping list. 
	* If the hopping already exists, add it as a different hopping with the same tag.
	*
	* @param[in] cellID Destintation cell ID.
	* @param[in] value  Hopping Kind.
	* @param[in] edge   Initial and final orbital ID.
	*/
	void AddHopping(const string& tag ,const hopping_t& hop);


	/**
	* Scale the hopping defined by its Destination cell ID and its edge. 
	* If it does not exists, nothing happens.
	*
	* @param[in] cellID Destintation cell ID.
	* @param[in] value  Value used for rescaling the hopping hop(cellID,edge).
	* @param[in] edge   Initial and final orbital ID.
	*/
	void ScaleHopping(cellID_t _cellID,value_t _value ,edge_t _edge);

	/**
	* Return the tag string from a given Destination Cell ID and edge. 
	*
	* @param[in] cellID Destintation cell ID.
	* @param[in] edge   Initial and final orbital ID.
	* @return Tag String defining the tag of the system.
	*/
	string get_tag(const hopping_list::cellID_t& cid,const hopping_list::edge_t edge)
	{
		string text_tag;
		for( auto& ti : cid)
			text_tag += to_string(ti)+" ";
		for( auto& ti : edge)
			text_tag += to_string(ti)+" ";
		return text_tag; 
	}

	/**
	* Parse the hoppings from a _hr.dat file
	* @param[in] WannierData Tuple containing the number of wannier basis, and strings read from the wannier_hr.dat file containing the wannier operator. 
	* @return hopping_list Returns a pointer to a modified hopping_list containing this hoppings.
	*/
	hopping_list& add_from_wannier( tuple<int, vector<string> > wannier_data  );


	/**
	* Add a set of hoppings with a random component from a hopping  file
	* @param[in] hopping_data Vector containingstrings read from the wannier_hr.dat file containing the wannier operator. 
	* @return hopping_list Returns a pointer to a modified hopping_list containing this hoppings.
	*/
	hopping_list& add_random_hopping_list( vector< vector<string> > disorder_data  );




	/**
	* Embedded the system into a supercell. This procedure changes the destination cell
	* of all hoppings, adding them when needed. It also changes the initial and final orbital ID in order to accomodate 
	* the inner atoms of the new supercell.
	* @param[in] cellDim New supercell over which the hoppings are defined. 
	* @return hopping_list Returns a pointer to a modified hopping_list defined over a new unit cells. 
	*/
	hopping_list wrap_in_supercell(const hopping_list::cellID_t& cellDim );


	/**
	* Transform the hopping list into CSR sparse matrix format, and save it to a file.
	* @param[in] output_filename Name of the file which will be used to save the hopping list . 
	*/
	void save_hopping_list_as_csr(string output_filename);


	//TO DO
    bool operator ==(hopping_list& y )
    {
		std::cout<<"EQUAL OPERATOR NOT IMPLEMENTED. RETURN TRUE"<<std::endl;
		return true;
	};
	
	private:
    int num_wann;
    cellID_t cellSizes;

    public:
    multimap< string  ,  hopping_t > hoppings;
    //map< string  ,  hopping_t > hoppings;
};







#endif
