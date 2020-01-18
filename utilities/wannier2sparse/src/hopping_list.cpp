#include "hopping_list.hpp"

hopping_list create_hopping_list( tuple<int, vector<string> > wannier_data  )
{
    hopping_list hl;
    hl.num_wann= get<0>(wannier_data); //number of wannier functions
    const vector<string> hopping_lines= get<1>(wannier_data); //strings containing the hopping data
    for (auto line : hopping_lines){
        stringstream ss(line);
        ss.precision( numeric_limits<double>::digits10+2);
       
        hopping_list::cellID_t cellID; 
        for( auto& x : cellID)
            ss>>x;

        hopping_list::edge_t vertex_edge; 
        for( auto& x : vertex_edge){ ss>>x; x-=1; }//The input is assumed to be zero based        
        
        double re,im; ss>>re>>im; 
        hopping_list::value_t hop_value(re,im);

        if( sqrt( norm(hop_value) )> numeric_limits<double>::epsilon() )
            hl.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop_value,vertex_edge) } );
    }
return hl; 
};   

hopping_list wrap_in_supercell(const hopping_list::cellID_t& cellDim,const hopping_list hl ){

    hopping_list sc_hl;;  //the hopping list for the supercell
    sc_hl.SetWannierBasisSize(hl.WannierBasisSize()); //increase the supercell dimension;
    sc_hl.SetBounds(cellDim); //increase the supercell dimension;

    //REPLICATE THESE HOPPINGS IN THE SUPERCELL
    hopping_list::cellID_t cellShift; 
    for(cellShift[2]=0; cellShift[2]< cellDim[2]; cellShift[2]++)
    for(cellShift[1]=0; cellShift[1]< cellDim[1]; cellShift[1]++)
    for(cellShift[0]=0; cellShift[0]< cellDim[0]; cellShift[0]++)
    {
        //Go through all the hopping list defined in the unit cell 
        //map the cell indexes into the new dimensions of the super cell
        //and add the values when it correspond
        for (auto const& key_hop : hl.hoppings){
            const auto hop = key_hop.second;

            auto cellID = get<0>(hop);
            auto value  = get<1>(hop);
            auto edge   = get<2>(hop);

             //Shift the cell given by the hopping vector
            //cellID = cellID + cellShift ;
            std::transform( cellID.begin(),cellID.end(),cellShift.begin(),cellID.begin(),  std::plus<int>() );

            //wrap tag_indexes around the super_cell 
            for( size_t i=0; i < cellID.size(); i++)
                cellID[i]=( cellID[i]+cellDim[i])%cellDim[i];

            //Shift both the beggining and end of the edge.
            edge[0] += index_aliasing(cellShift,cellDim)*hl.WannierBasisSize(); //This is assume to be at the origin and shifted by cellShift
            edge[1] += index_aliasing(cellID,cellDim)*hl.WannierBasisSize(); //This is originally shidted by cellID and one added an aditional shift by cellShift.
            cellID = {0,0,0}; //After wrapping, all atoms belong to the supercell, therefore, cellID is always zero.

            //Insert the hoppings
            assert( edge[0]< sc_hl.WannierBasisSize() );
            assert( edge[1]< sc_hl.WannierBasisSize() );
            assert( cellID== hopping_list::cellID_t({0,0,0}) );
            string cellID_tag = get_tag(cellID,edge);

            if (sc_hl.hoppings.count(cellID_tag) ==  0)
                sc_hl.hoppings.insert( { cellID_tag, hopping_list::hopping_t(cellID,value,edge) } );
            else
                get<1>(sc_hl.hoppings[cellID_tag]) += value; 
        }
    }
    sc_hl.SetBounds(hopping_list::cellID_t({1,1,1})); //The wrapped cell is bounded in itself
return sc_hl; 
};

void hopping_list::AddHopping(hopping_list::cellID_t _cellID,hopping_list::value_t _value ,hopping_list::edge_t _edge)
{
	auto tag = get_tag(_cellID,_edge);
	if( this->hoppings.count(tag) > 0 )
	{
		get<1>( this->hoppings[ tag ] )+=_value; //this add sval*value to the spincurrent value list;
		return ;
	}
	this->hoppings.insert( {tag, hopping_t(_cellID,_value,_edge) } );
	return ;
};




void save_hopping_list_as_csr(string output_filename,const hopping_list& hl)
{
    const size_t dim = hl.WannierBasisSize();
    SparseMatrix_t output(dim,dim);
    std::vector<Triplet_t> coefficients;            // list of non-zeros coefficients
	const double CUTOFF = 0;
	for(auto const& elem : hl.hoppings)
    {
        const auto value  = get<1>(elem.second);
        const auto edge   = get<2>(elem.second);
        if( (value.real() != 0 || value.imag()!= 0)  && edge[1]>= edge[0] )
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

