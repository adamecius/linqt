#include "hopping_list.hpp"

//Temporal time measuring
  

void hopping_list::AddHopping(	const cellID_t& _cellID, const value_t& _value , const edge_t& _edge)
{
	const auto tag = get_tag(_cellID,_edge);
	const hopping_kind hop_kind(_value);
	const hopping_t hopping(_cellID,hop_kind,_edge);
	this->hoppings.insert( {tag,hopping} );
	return ;
};


void hopping_list::AddHopping(	const cellID_t& _cellID,const hopping_kind& _hop , const edge_t& _edge)
{
	const auto tag = get_tag(_cellID,_edge);
	const hopping_t hopping(_cellID,_hop,_edge);
	this->hoppings.insert( {tag,hopping} );
	return ;
};

void hopping_list::AddHopping(const string& tag ,const hopping_t& hop)
{
	this->hoppings.insert( {tag, hop} );
	return ;
};

hopping_list& hopping_list::add_random_hopping_list( vector< vector<string> > disorder_data  )
{

	for( auto group : disorder_data )
	{

		//THE LIST IS PRESENTED IN REVERSED ORDER
		//INPUT SHOULD BE:
		//distributiton type
		// x x x	 //shift kind 
		// xio xif yio yif  //orbital pair x and y
		// Ret_min Imt_min Ret_max Imt_max concentration
		
		//READ CONCENTRATION
		double concentration;
		{
			auto line = group.back(); std::cout<<"concentration: "<<line<<std::endl;
			concentration =  std::stod(line);
			group.pop_back();		
		}
		
		//READ MAX_MIN_VALUES
		hopping_list::value_t min_val,max_val;
		{
			auto line = group.back();
			string::size_type sz;     // alias of size_t
			min_val.real( stod(line,&sz) ); line = line.substr(sz);
			min_val.imag( stod(line,&sz) ); line = line.substr(sz);
			max_val.real( stod(line,&sz) ); line = line.substr(sz);
			max_val.imag( stod(line,&sz) ); line = line.substr(sz);
			std::cout<<"min,max: "<<min_val<<","<<max_val<<std::endl;
			assert( std::norm(min_val)<=std::norm(max_val) );
			group.pop_back();			
		}
  
		//READ ORBITAL_LIST
		vector< hopping_list::edge_t > vertexes;
		{
			auto line = group.back();
			string::size_type sz;     // alias of size_t
			int i0,i1;
			while( line.size()!= 0 ) 
			{
				i0 = std::stoi(line,&sz); line = line.substr(sz); 
				i1 = std::stoi(line,&sz); line = line.substr(sz); std::cout<<"i0,i1: "<<i0<<","<<i1<<std::endl;
				vertexes.push_back( hopping_list::edge_t({i0-1,i1-1}) ); //Assume 1 based
				
			}
			group.pop_back();		
		}

		//READ CELL_ID
		hopping_list::cellID_t cellID;
		{
			auto line = group.back();
			string::size_type sz;     // alias of size_t
			for( auto& x : cellID)
			{
				x = std::stoi(line,&sz); line = line.substr(sz); 
			}
			group.pop_back();		
		}
		for( auto line : group )
			std::cout<<line<<std::endl;
		//READ DISTRIBUTION
		string distribution;
		{
			auto line = group.back();
			distribution = line;
			group.pop_back();		
		}
		//Create a hopping which is common to all the orbitals
		auto hopping = hopping_kind().Concentration(concentration)
									 .DistributionType(distribution,min_val,max_val);
	
		for( auto vertex : vertexes )
			this->AddHopping( cellID,hopping,vertex );
	}
 
return *this; 
};   



hopping_list& hopping_list::add_from_wannier( tuple<int, vector<string> > wannier_data  )
{
    this->SetWannierBasisSize( get<0>(wannier_data) ); //number of wannier functions
    const vector<string> hopping_lines= get<1>(wannier_data); //strings containing the hopping data
    
    std::vector< hopping_t > hoppings;
	hoppings.reserve( hopping_lines.size() );
	double avg_hop_norm = 0;
    for (auto line : hopping_lines){
        stringstream ss(line);
        ss.precision( numeric_limits<double>::digits10+2);
       
        hopping_list::cellID_t cellID; 
        for( auto& x : cellID)
            ss>>x;

        hopping_list::edge_t vertex_edge; 
        for( auto& x : vertex_edge){ ss>>x; x-=1; }//The input is assumed to be zero based        
        
        double re,im; ss>>re>>im;
        hopping_list::value_t hop(re,im) ;
        avg_hop_norm += sqrt(std::norm(hop));
        hopping_kind value(hop); 
		hoppings.push_back( hopping_t(cellID,value,vertex_edge) );	
    }
    avg_hop_norm /= hopping_lines.size();
	assert(avg_hop_norm>0);

    for (auto elem : hoppings )
    {
		auto  cellID = get<0>(elem);
		auto  value  = get<1>(elem);
		auto  edge   = get<2>(elem);

		auto hop_norm=sqrt(std::norm(value()));
    	if( hop_norm/avg_hop_norm > 1e-10 )  
			this->AddHopping(cellID,value,edge);
	}
    

return *this; 
};   


hopping_list hopping_list::wrap_in_supercell(const hopping_list::cellID_t& cellDim ){

    hopping_list sc_hl;  //the hopping list for the supercell
    const auto WBB = this->WannierBasisSize();
    sc_hl.SetWannierBasisSize(WBB); //increase the supercell dimension;
    sc_hl.SetBounds(cellDim); //increase the supercell dimension;

	cout<<"\nINVOCATION OF FUNCTION: wrap_in_supercell. "<<endl
		<<"The cell size is ("<<cellDim[0]<<","<<cellDim[1]<<","<<cellDim[2]<<")"<<endl
		<<"In this particular instance, a minimum of "
		<<sizeof(hopping_t)* this->hoppings.size()* cellDim[0]* cellDim[1]* cellDim[2]*1e-9
		<<" GB is required.\n"<<endl
		<<"The insertion of a hopping in the supercell ";


    //REPLICATE THESE HOPPINGS IN THE SUPERCELL
    hopping_list::cellID_t cellShift; 
	//Go through all the hopping list defined in the unit cell 
	//map the cell indexes into the new dimensions of the super cell
	//and add the values when it correspond
	auto start = high_resolution_clock::now(); 
	double hopping_insertion_time = 0.0;
	double num_hops = this->hoppings.size();
	for (auto const& key_hop : this->hoppings)
	{
		const auto _tag		= key_hop.first;
		const auto _hop		= key_hop.second;
		const auto _cellID	= get<0>(_hop);
		const auto _value	= get<1>(_hop);
		const auto _edge	= get<2>(_hop);

		auto value	= _value;
	    for(cellShift[2]=0; cellShift[2]< cellDim[2]; cellShift[2]++)
	    for(cellShift[1]=0; cellShift[1]< cellDim[1]; cellShift[1]++)
	    for(cellShift[0]=0; cellShift[0]< cellDim[0]; cellShift[0]++)
	    {
			auto cellID	= _cellID; 
			auto edge	= _edge;
	
			//Shift the cell given by the hopping vector
	        //cellID = cellID + cellShift ;
            //wrap tag_indexes around the super_cell 
            for( size_t i=0; i < cellID.size(); i++)
                cellID[i]=( cellID[i]+cellShift[i]+cellDim[i])%cellDim[i];

	            //Shift both the beggining and end of the edge.
	            edge[0] += index_aliasing(cellShift,cellDim)*WBB; //This is assume to be at the origin and shifted by cellShift
	            edge[1] += index_aliasing(cellID,cellDim)*WBB; //This is originally shidted by cellID and one added an aditional shift by cellShift.
	            cellID = hopping_list::cellID_t({0,0,0}) ; //After wrapping, all atoms belong to the supercell, therefore, cellID is always zero.
	
	            //Insert the hoppings
	            assert( edge[0]< sc_hl.WannierBasisSize() );
	            assert( edge[1]< sc_hl.WannierBasisSize() );         
	            sc_hl.AddHopping(cellID, value() , edge ) ;
		}
		if( hopping_insertion_time == 0.0)
		{
			auto stop = high_resolution_clock::now(); 
			auto duration = duration_cast<microseconds>(stop - start); 
			hopping_insertion_time =  duration.count()*1e-6;
			cout <<"takes "<< hopping_insertion_time <<" seconds."
				 <<" The stimated completion time: "<<  hopping_insertion_time*num_hops <<" seconds."<< endl;
		}
    }
    sc_hl.SetBounds(hopping_list::cellID_t({1,1,1})); //The wrapped cell is bounded in itself
return sc_hl; 
};


void hopping_list::save_hopping_list_as_csr(string output_filename)
{
    const size_t dim = this->WannierBasisSize();
    SparseMatrix_t output(dim,dim);
    std::vector<Triplet_t> coefficients;            // list of non-zeros coefficients
	const double CUTOFF = 0;

	//Compute the average hoppinh value
	double avg_hop_norm=0.0;
	for(auto const& elem : this->hoppings)
    {
		auto hop   = get<1>(elem.second);
		auto value = sqrt(std::norm(hop()));
		avg_hop_norm += value; 
    }
	avg_hop_norm /= this->hoppings.size();


	for(auto const& elem : this->hoppings)
    {

		auto hop   = get<1>(elem.second);
        auto edge  = get<2>(elem.second);
		auto value = hop();
		auto abs_norm = sqrt(std::norm(value));

        if( abs_norm/avg_hop_norm > 1E-10   && edge[1]>= edge[0] )
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

