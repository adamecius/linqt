#include "lattice.hpp"

namespace NumCal {

void Lattice::PrintIrregularOnsiteProfile(std::string output )
{
	std::ofstream output_file(output.c_str());

	for ( integer k=0; k< IrregHam().Matrix().outerSize(); ++k)
		for (  my::SparseMatrix::InnerIterator it(IrregHam().Matrix() ,k); it; ++it)
			if( it.row()== it.col() )
			{
				integer i0,i1,i2,io,is;
				IndexToIndexes(it.col(),i0,i1,i2,io,is);
				const  real r[2]={
						A[0][0]*i0+A[1][0]*i1 + io*Delta[0],
						A[0][1]*i0+A[1][1]*i1 + io*Delta[1]
				};
				if(is==0)
					output_file<<r[0]<<" "<<r[1]<<" "<<it.value().real()<<NumCal::endl;
			}
	output_file.close();

}

void Lattice::RescaleHamiltonian(const real _Emin, const real _Emax,
		const real _cutoff) {
	if (IsHamiltonianSetted()) {
		real Energy_scal = (_Emax - _Emin) / (2. * _cutoff);
		real Energy_shift = (_Emax + _Emin) / 2.;
		hop_list_.ShiftOnSite(Energy_shift);
		hop_list_.RescaleHopping(1 / Energy_scal);
		IrregHam().Rescale(_Emin, _Emax, _cutoff);
	} else {
		NumCal::cerr << "Error while rescaling the hamiltonian,"
				<< "the hamiltonian is not setted" << NumCal::endl;
		std::exit(-1);
	}
}

real Lattice::SetVolume() {
	volume_ = -LatticeVector(0)[2] * LatticeVector(1)[1] * LatticeVector(2)[0]
			+ LatticeVector(0)[1] * LatticeVector(1)[2] * LatticeVector(2)[0]
			+ LatticeVector(0)[2] * LatticeVector(1)[0] * LatticeVector(2)[1]
			- LatticeVector(0)[0] * LatticeVector(1)[2] * LatticeVector(2)[1]
			- LatticeVector(0)[1] * LatticeVector(1)[0] * LatticeVector(2)[2]
			+ LatticeVector(0)[0] * LatticeVector(1)[1] * LatticeVector(2)[2];

	if (volume_ == 0) {
		NumCal::cerr
				<< "The volume cannot be zero, if you are trying to perform"
				<< "a 2D or 1D simulation please consider setting the unused dimension vector to the unitary "
				<< NumCal::endl;
		std::exit(-1);
	}

	return volume_;
}

void Lattice::SetReciprocalLatticeVectors() {
	SetVolume();

	rec_lat_[0][0] = -(LatticeVector(1)[2] * LatticeVector(2)[1]
			- LatticeVector(1)[1] * LatticeVector(2)[2]) * 2. * M_PI / Volume();
	rec_lat_[0][1] = (LatticeVector(1)[2] * LatticeVector(2)[0]
			- LatticeVector(1)[0] * LatticeVector(2)[2]) * 2. * M_PI / Volume();
	rec_lat_[0][2] = -(LatticeVector(1)[1] * LatticeVector(2)[0]
			- LatticeVector(1)[0] * LatticeVector(2)[1]) * 2. * M_PI / Volume();

	rec_lat_[1][0] = (LatticeVector(0)[2] * LatticeVector(2)[1]
			- LatticeVector(0)[1] * LatticeVector(2)[2]) * 2. * M_PI / Volume();
	rec_lat_[1][1] = -(LatticeVector(0)[2] * LatticeVector(2)[0]
			- LatticeVector(0)[0] * LatticeVector(2)[2]) * 2. * M_PI / Volume();
	rec_lat_[1][2] = (LatticeVector(0)[1] * LatticeVector(2)[0]
			- LatticeVector(0)[0] * LatticeVector(2)[1]) * 2. * M_PI / Volume();

	rec_lat_[2][0] = -(LatticeVector(0)[2] * LatticeVector(1)[1]
			- LatticeVector(0)[1] * LatticeVector(1)[2]) * 2. * M_PI / Volume();
	rec_lat_[2][1] = (LatticeVector(0)[2] * LatticeVector(1)[0]
			- LatticeVector(0)[0] * LatticeVector(1)[2]) * 2. * M_PI / Volume();
	rec_lat_[2][2] = -(LatticeVector(0)[1] * LatticeVector(1)[0]
			- LatticeVector(0)[0] * LatticeVector(1)[1]) * 2. * M_PI / Volume();

}

void Lattice::IndexToIndexes(const my::integer k0, my::integer& i0,
		my::integer& i1, my::integer& i2, my::integer& o0, my::integer& s0) {
	s0 = k0 % (SpinNumber());
	o0 = (k0 / SpinNumber()) % (OrbitalNumber());
	i2 = (k0 / SpinNumber() / OrbitalNumber()) % (CellsInDir(2));
	i1 = (k0 / SpinNumber() / OrbitalNumber() / CellsInDir(2))
			% (CellsInDir(1));
	i0 = (k0 / SpinNumber() / OrbitalNumber() / CellsInDir(2) / CellsInDir(1))
			% (CellsInDir(0));
}
;

void Lattice::ReadHoppingList() {

	NumCal::cout << "-------------READING FROM HOPPING FILES-------------"
			<< NumCal::endl << NumCal::endl;

	std::string hop_filename(label_ + ".hop");
	std::ifstream hop_file(hop_filename.c_str());

	if (!hop_file.is_open()) {
		NumCal::cerr << "The hopping file:" << hop_filename
				<< " was not found, please submit one" << NumCal::endl;
		std::exit(-1);
	} else
		NumCal::cout << "Reading from " << hop_filename << " :" << NumCal::endl;

	std::vector<integer> Dr(3);
	integer o0, o1, s0, s1;
	real reVal, imVal;

	//Read the total number of lines in the file
	integer line_num;
	{
		std::string line;
		for (line_num = 0; std::getline(hop_file, line); ++line_num)
			;
	}

	//Go to the begining of the file
	hop_file.clear();
	hop_file.seekg(0, std::ios::beg);

	//Reading the hoppings
	for (integer l = 0; l < line_num; l++) {
		hop_file >> o0 >> s0 >> Dr[0] >> Dr[1] >> Dr[2] >> o1 >> s1 >> reVal
				>> imVal;
		integer init_idx = InnerCellIndex(o0, s0);
		integer final_idx = InnerCellIndex(o1, s1);
		scalar val(reVal, imVal);
		hop_list_.AddHopping(init_idx, Hopping(final_idx, o1, s1, Dr, val));
	}
	//Print the Hopping list
	hop_list_.PrintHoppingList();

}

}
