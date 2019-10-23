/*
 * sparse_matrix.hpp
 *
 *  Created on: Aug 7, 2016
 *      Author: jgarcia
 */

#ifndef SPARSE_MATRIX_HPP_
#define SPARSE_MATRIX_HPP_

//***********Standard C/C++ functions****************/
///Defines the std::vector class
#include <vector>
///Defines the Eigen::SparseMatrix and Eigen::Triplet  classes
#include "Eigen/Sparse"
/// Define the different types of the program
#include "types_definitions.hpp"

///Namespace for the alias given to some of the objetcs of the EigenClass
namespace sparse
{
typedef Eigen::SparseMatrix< my::scalar, Eigen::RowMajor, my::integer > EigenMat;
typedef Eigen::Triplet<my::scalar> Entry;
}

///Namespace for the classes created for this project
namespace my {
  ///The Sparse Matrix Class. Inheritated from sparse::EigenMat
  class SparseMatrix : public sparse::EigenMat
  {
  public:
    //******************CONSTRUCTORS***************************/
    ///The Null constructor
    SparseMatrix(): sparse::EigenMat(), matrix_is_set_(false){}
    ///The Default constructor
    SparseMatrix( const my::integer _ncol,
		  const my::integer _nrow):
		    sparse::EigenMat(_ncol,_nrow),
		    matrix_is_set_(false)
    {}

    //***************PUBLIC METHODS***************************/
    /// Reserve an stimated ammount of memory for the sparse matrix
    void
    Reserve(const my::integer _size)
    {
      matEntry_.reserve(_size);
    }

    ///Returns the total number of raw entries in the matrix
    /*! Returns the total number of raw entries in the matrix.
     *  This method will count repeated and non-zero entries
     */
    my::integer
    NumRawEntries() const
    {
      return matEntry_.size();
    }

    ///Add a new Raw entry in the matrix
    /*! This method will append a new entry to the entry list,
     * not making any check on the previous entries
     */
    void
    AddNewRawEntry(const sparse::Entry _triplet)
    {
      matEntry_.push_back(_triplet);
    }

    ///Set the sparse matrix using the raw_entries.
    /*! The list of entries after this operation is destroyed
     *
     */
    void
    SetFromRawEntries()
    {
      this->setFromTriplets(matEntry_.begin(), matEntry_.end() );
      matrix_is_set_=true;
      matEntry_= std::vector<sparse::Entry>(0);
  	for(int i=0;i<matEntry_.size();i++ )
  		std::cout<<matEntry_[i].col()<<" "<<matEntry_[i].row()<<" "<<matEntry_[i].value()<<" "<<std::endl;

    }

    //*************GETTERS AND SETTERS***************************/
    sparse::Entry
    RawEntry(const my::integer _idx) const
    {
      return matEntry_[_idx];
    }

    sparse::Entry&
    RawEntry(const my::integer _idx)
    {
      return matEntry_[_idx];
    }

    bool
    IsSet() const
    {
      return matrix_is_set_;
    }
private:
	std::vector<sparse::Entry> matEntry_;
	bool matrix_is_set_;
};

};
#endif /* SPARSE_MATRIX_HPP_ */
