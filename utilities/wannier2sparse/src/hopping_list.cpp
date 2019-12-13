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

            //Insert the hoppinhs
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
