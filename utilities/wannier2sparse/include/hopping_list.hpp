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
#include<iostream>
#include<limits>
#include<algorithm>
#include"sparse_matrix.hpp"
#include <iostream>

using namespace std;

inline int 
index_aliasing(const array<int, 3>& index,const array<int, 3>& bound )
{
    return ( (index[2]+bound[2])%bound[2] * bound[1] + (index[1]+bound[1])%bound[1] ) * bound[0] + (index[0]+bound[0])%bound[0] ;
}


struct hopping_list
{
    typedef complex<double> value_t;
    typedef array<int, 2> edge_t;
    typedef array<int, 3> cellID_t;
    typedef tuple< cellID_t,value_t,edge_t > hopping_t;

    hopping_list():cellSizes({1,1,1}), num_wann(0){};

    inline int WannierBasisSize() const
    {
        return this-> num_wann ; 
    };

    inline void SetWannierBasisSize(const int num_wann)
    {
        assert( num_wann > 0 );
        this-> num_wann  = num_wann;
        return ; 
    };

    inline void SetBounds(const cellID_t& cellSizes)
    {
        assert( this-> num_wann > 0 && cellSizes[0]>0&& cellSizes[1]>0&&cellSizes[2]>0 );
        this->cellSizes = cellSizes;
        for(const auto& x: cellSizes ) 
            this->num_wann *=x;
        return ; 
    };

    cellID_t Bounds()
    {
        return array<int, 3>({cellSizes[0],cellSizes[1],cellSizes[2]});
    };

    int cellID_index(const cellID_t cidx )
    {
        const cellID_t bounds = this->Bounds();
       
        int index = 0;
        for( int i = 0; i+1 < cidx.size(); i++ )
            index += ( (cidx[i]+bounds[i])%bounds[i] )*bounds[i+1];
        index += ( cidx.back()+bounds.back() )%( bounds.back() );
        return index;
    };
    bool operator ==(hopping_list& y )
    {
        bool list_equal = true;
        for( auto const& elem: y.hoppings )
        {
            auto key = elem.first;
            if( this->hoppings.count(key) == 0 )
            {
                std::cout<<"The key: "<<key<<" was not found when comparing the hopping lists"<<std::endl;
                list_equal = false;    
                break;
            }
            list_equal*= (bool)(get<0>(this->hoppings[key])==get<0>(elem.second));
            list_equal*= (bool)(get<2>(this->hoppings[key])==get<2>(elem.second));

            auto val_diff= (get<1>(this->hoppings[key])-get<1>(elem.second ))/2.0;
            auto val_sum = (get<1>(this->hoppings[key])+get<1>(elem.second ))/2.0;

            if(!( val_diff.real()==0&& val_diff.imag()==0 ) )
            {
                list_equal*= (bool)(
                                std::fabs(val_diff.real()/val_sum.real()) < std::numeric_limits<double>::epsilon() && 
                                std::fabs(val_diff.imag()/val_sum.imag()) < std::numeric_limits<double>::epsilon() 
                                );
                if(!list_equal)
                    std::cout<<"The keys "<<key<<" have hoppings with a percentile difference higher than "<< std::numeric_limits<double>::epsilon()<<std::endl;
            }
        }
        return  ( this->num_wann ==y.num_wann)&&
                (this->cellSizes==y.cellSizes)&&
                list_equal;
    }

    int num_wann;
    cellID_t cellSizes;
    map<string, hopping_t > hoppings; 
};

hopping_list create_hopping_list( tuple<int, vector<string> > wannier_data  );

hopping_list wrap_in_supercell(const hopping_list::cellID_t& cellDim,const hopping_list hl);

inline string get_tag(const hopping_list::cellID_t& cid,const hopping_list::edge_t edge){
    string text_tag;
    for( auto& ti : cid)
        text_tag += to_string(ti)+" ";
    for( auto& ti : edge)
        text_tag += to_string(ti)+" ";
return text_tag; 
}

inline array<int,5> tag_to_indices(const string& tag){
    array<int,5> indices;
    stringstream ss(tag);
    for( auto& ti : indices)
        ss>>ti;
return indices; 
}

inline void save_hopping_list_as_csr(string output_filename,const hopping_list& hl)
{
    const size_t dim = hl.WannierBasisSize();
    SparseMatrix_t output(dim,dim);
    std::vector<Triplet_t> coefficients;            // list of non-zeros coefficients
    for(auto const& elem : hl.hoppings)
    {
        const auto value  = get<1>(elem.second);
        const auto edge   = get<2>(elem.second);
        coefficients.push_back(Triplet_t(edge[0],edge[1],value) );
    }
    output.setFromTriplets(coefficients.begin(), coefficients.end());
    output.makeCompressed();


	std::ofstream matrix_file ( output_filename.c_str()) ;

	//READ DIMENSION OF THE MATRIX
	matrix_file<<dim<<" "<<output.nonZeros()<<std::endl; 

    //save values first
    for (int k=0; k<output.outerSize(); ++k)
    for (SparseMatrix_t::InnerIterator it(output,k); it; ++it)
        matrix_file<<it.value().real()<<" "<<it.value().imag()<<" ";
    matrix_file<<std::endl;

    //save the columns
    for (int k=0; k<output.outerSize(); ++k)
    for (SparseMatrix_t::InnerIterator it(output,k); it; ++it)
        matrix_file<<it.index()<<" ";
    matrix_file<<std::endl;

    //save the indices to columns
    for (int k=0; k<output.outerSize()+1; ++k)
        matrix_file<<*( output.outerIndexPtr() + k ) <<" ";
    matrix_file<<std::endl;

    matrix_file.close();

return ; 
}

#endif
