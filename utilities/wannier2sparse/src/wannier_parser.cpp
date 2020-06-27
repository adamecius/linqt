#include "wannier_parser.hpp"


vector< vector<string> > read_disorder_file(const string disorder_filename)
{
	ifstream input_file(disorder_filename.c_str()); 
    input_file.precision( numeric_limits<double>::digits10+2);
    if(! input_file.is_open())
		return vector< vector<string> >();
    else
        std::cout<<"Found disorder file:"<<disorder_filename<<std::endl;

	//Get all lines into the memory
    string line; 
	vector<string>  lines;
    for( int counter = 0; getline(input_file, line); counter++)
		lines.push_back( line );


	//Group lines depending on if there is a blank space between them
	vector<string> group;
	vector< vector<string> >  line_groups;
    for( auto line : lines  )
    {
		if( line.size() == 0) 
		{
			if( group.size() != 0 )
			{
				line_groups.push_back(group);
				group = vector<string>();
			}
		}
		else
			group.push_back( line );
	}
	if( group.size() != 0 )
		line_groups.push_back(group);


    input_file.close();        
	return line_groups;
};


tuple<int, vector<string> > read_wannier_file(const string wannier_filename)
{
    int num_wann = 0;   // the number of Wannier functions
    int nrpts = 0;      //number of Wigner-Seitz grid-points
    vector< string > wz_gpoints; 
    vector< string > hopping_list; 

    ifstream input_file(wannier_filename.c_str()); 
    assert(input_file.is_open());
    input_file.precision( numeric_limits<double>::digits10+2);

    string line; 
    for( int counter = 0; getline(input_file, line); counter++){
        switch( counter ){
            case 0 : continue;  //ignore data or comments
            
            case 1 : num_wann = safe_stoi(line  ); assert(num_wann>0); continue;
            
            case 2 : nrpts = safe_stoi(line); assert(nrpts>0); continue;
            
            case 3 :{   //Read all Wigner-Seitz grid-points
                const int nrpts_lines = nrpts/15;   //the format impose 15 grid-points per line.  
                wz_gpoints.push_back(line);         //Note the first line is always read so we are always reading nrpts_lines+1
                for(int n=0; n < nrpts_lines ; n++)
                {
                    getline(input_file, line);
                    wz_gpoints.push_back(line);
                    counter++;
                }
               continue;
            }
            default:
                hopping_list.push_back( line );
        }
    }
    input_file.close();        


return tuple<int, vector<string> >(num_wann,hopping_list); 
};

vector< tuple<string, array<double, 3> > > read_xyz_file(const string xyz_filename)
{
    typedef tuple<string, array<double, 3> > xyz_elem;
    vector< xyz_elem > xyz_data;
    std::cout<<"Opening File: "<<xyz_filename.c_str()<<std::endl;
    ifstream input_file(xyz_filename.c_str()); assert(input_file.is_open());
    input_file.precision( numeric_limits<double>::digits10+2);
    int num_sites;
    input_file>>num_sites; assert(num_sites>0);

    std::string label;
    array<double, 3>  pos;
    for( int i = 0; i < num_sites; i++)
    {
       input_file>>label>>pos[0]>>pos[1]>>pos[2];
       xyz_data.push_back(xyz_elem(label,pos) );
    } 
    return xyz_data;
};

array< array<double,3> , 3 >  read_unit_cell_file(const string uc_filename)
{
    constexpr int DIM = 3;
    ifstream input_file(uc_filename.c_str()); assert(input_file.is_open());
    input_file.precision( numeric_limits<double>::digits10+2);
    array< array<double,DIM> , DIM > unit_cell;
    for( auto & lat_vec : unit_cell )
    for( auto & li : lat_vec )
       input_file>>li;
    return unit_cell;

};


