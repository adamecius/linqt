#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <string>
#include <vector> // used for std::vector
#include <fstream> // used for std::ifstream and std::ofstream
#include <iostream>//used for NumCal::cout and std::err
//Defines the types in the program
#include "types_definitions.hpp"
//Defines the irregular hamiltonian class
#include "irregular_hamiltonian.hpp"
//Defines the lattice geometry, in this case honeycomb
#include "lattice_geometry.h"
//Defines the hopping class
#include "hopping.hpp"
//Defines the hopping class
#include "hopping_list.hpp"
//Defines the hopping class
#include "efficient_mod.hpp"
//Defines the fourier transform functions class
#include "fourier_transform.hpp"

#include "kpm_parallel.hpp"

#include "mpi_util.hpp"

namespace NumCal {
/** \addtogroup lattice
 *  @{
 */

/**\class Lattice
 *
 * \brief A multi-purpose tight-binding lattice class
 *
 * This class is aimed to served as a multipurpose tight-binding
 * lattice class, which can be used to store different model for
 * tight-binding Hamiltonians of solid-state systems. This class
 * read the hoppings from a .hop file and construct a regular hopping
 *
 *  It st geometrical properties, hoppings
 * and other variables. It serve as  a host of the irregular and regular
 * part of the hamiltonian also.
 *
 * \note This class is still in development process, use it at your own risk!
 *
 * \version 1.0
 *
 */
class Lattice {
public:
	/********************************CONSTRUCTORS*****************************/
	///The Main constructor
	/*!
	 * This constructor will try to open the file _config_filename
	 * and will read the following parameters:\n
	 * The label given to the simulation.\n
	 * The lattice vectors: lat0, lat1, lat2.\n
	 * The number of unit cells in lat_0 direction\n
	 * The number of unit cells in lat_1 direction\n
	 * The number of unit cells in lat_2 direction\n
	 * The number of orbital per unit cell.\n
	 * The total spin number of the simulation
	 * The coordination number.\n
	 */

	Lattice(const std::string _config_filename) :
			hamiltonian_setted_(false) {

		//Initialize the lattice vectors
		lat_ = std::vector < std::vector<real> > (3);
		for (short i = 0; i < 3; ++i)
			lat_[i] = std::vector < real > (3, 0);

		//The CellInDir Vector
		nsite_ = std::vector < integer > (3);

		//Initialize the reciprocal lattice vectors
		rec_lat_ = std::vector < std::vector<real> > (3);
		for (short i = 0; i < 3; ++i)
			rec_lat_[i] = std::vector < real > (3, 0);

		//Opens the config file
		std::ifstream config_file(_config_filename.c_str());

		///Look for the line lattice_info
		std::string line;
		bool found = false;
		for (int line_num = 0; std::getline(config_file, line); ++line_num)
			if (line == "lattice_info") {
				found = true;
				break;
			}

		if (found) {
			config_file >> label_ >> lat_[0][0] >> lat_[0][1] >> lat_[0][2]
					>> lat_[1][0] >> lat_[1][1] >> lat_[1][2] >> lat_[2][0]
					>> lat_[2][1] >> lat_[2][2] >> nsite_[0] >> nsite_[1]
					>> nsite_[2] >> norb_ >> spin_ >> coord_;
			config_file.close();

			SetReciprocalLatticeVectors();

		} else {
			config_file.close();
			NumCal::cerr
					<< "The config file does not posses a lattice_info section. The simulation cannot proceed"
					<< NumCal::endl;
			std::exit(-1);
		}
		Hirr_ = IrregularHamiltonian(TotalOfOrbitals(), TotalOfOrbitals());
	}

	/********************************GETTERS*****************************/

	///Returns the Lattice Label
	inline std::string Label() const {
		return label_;
	}
	///Returns the Total number of sites
	inline integer TotalOfOrbitals() const {
		return nsite_[0] * nsite_[1] * nsite_[2] * SpinNumber()
				* OrbitalNumber();
	}
	;

	///Returns the ith-lattice vector
	inline std::vector<real> LatticeVector(const integer i) const {
		return lat_[i];
	}
	;

	///Returns the ith-reciprocal lattice vector
	inline std::vector<real> ReciprocalLatticeVector(const integer i) const {
		return rec_lat_[i];
	}
	;

	///Returns the total number of cells in the direction of the ith-lattice vector
	inline integer CellsInDir(const integer i) const {
		return nsite_[i];
	}
	;

	///Returns the total number of orbitals per unit cell
	inline integer OrbitalNumber() const {
		return norb_;
	}
	;

	///Returns the Spin number
	inline integer SpinNumber() const {
		return spin_;
	}
	;

	///Returns the coordination number
	inline integer CoordinationNumber() const {
		return coord_;
	}
	;

	///Returns a boolean indicating if the hamiltonian is setted or not
	bool IsHamiltonianSetted() const {
		return hamiltonian_setted_;
	}

	///Returns the irregular part of the Hamiltonian
	IrregularHamiltonian&
	IrregHam() {
		return Hirr_;
	}

	///returns the volume of the sample
	real Volume() const {
		return volume_;
	}

	/********************************SETTERS*****************************/

	///Set the volume of the sample
	real SetVolume();

	///Set the reciprocal lattice vectors
	void SetReciprocalLatticeVectors();

	/********************************Index Conversion Utilities*****************************/

	///Convert the lattice indexes into the site index
	inline integer IndexesToIndex(const integer i0, const integer i1,
			const integer i2, const integer o0, const integer s0) {
		return CellIndex(i0, i1, i2) + InnerCellIndex(o0, s0);
	}

	///Convert the site index into the lattice indexes
	void
	IndexToIndexes(const integer k0, integer& i0, integer& i1, integer& i2,
			integer& o0, integer& s0);

	inline integer CellIndex(const integer i0, const integer i1,
			const integer i2) {
		return (( //nat2
		( // nat 1
		EffMod(i0, CellsInDir(0)) //nat0
		) * CellsInDir(1) + EffMod(i1, CellsInDir(1))) * CellsInDir(2)
				+ EffMod(i2, CellsInDir(2))) * OrbitalNumber() * SpinNumber();
	}

	inline integer InnerCellIndex(const integer o0, const integer s0) {
		return o0 * SpinNumber() + s0;
	}

	///Reserve a estimated memory for the irregular part of the hamiltonian
	void ReserveIrrHamSpace(const integer _stimated_entries) {
		IrregHam().Reserve(_stimated_entries);
	}

	///Initialize the final version of the hamiltonian in the memory.
	bool SetTotalHamiltonian() {
		IrregHam().Refine();
		hamiltonian_setted_ = true;
		return hamiltonian_setted_;
	}

	///Rescale the Hamiltonian between (-cutoff, cutoff)
	/*!
	 * This function is of special interest for the KPM method. It
	 * is important to notice that it rescale both the regular
	 * and irregular hamiltonian and both have to be setted.
	 */
	void RescaleHamiltonian(const real _Emin, const real _Emax,
			const real _cutoff);

	/// Priknt disorder onsite
	void PrintIrregularOnsiteProfile(std::string output);

	///Read the list of hoppings from the hopping file
	void ReadHoppingList();

	void ApplyHamiltonian(const integer memSep, const scalar alpha,
			const scalar* x, const scalar beta, scalar* y) {

	//	#pragma omp parallel  for
		for (integer i = 0; i < memSep * TotalOfOrbitals(); i += memSep)
			y[i] = beta * y[i];

	//#pragma omp parallel for 	//This is only for small numbers of cores
		for (int i0 = 0; i0 < CellsInDir(0); i0++)
			for (int i1 = 0; i1 < CellsInDir(1); i1++)
				for (int i2 = 0; i2 < CellsInDir(2); i2++)
					for (int iorbs = 0; iorbs < hop_list_.TotNumRegHops();
							iorbs++)
						for (integer h = 0;
								h < hop_list_.HopOfOrbs(iorbs).size(); h++) {
							const integer j0 = i0
									+ hop_list_.HopOfOrbs(iorbs)[h].Dr[0];
							const integer j1 = i1
									+ hop_list_.HopOfOrbs(iorbs)[h].Dr[1];
							const integer j2 = i2
									+ hop_list_.HopOfOrbs(iorbs)[h].Dr[2];

							const integer row = memSep
									* (CellIndex(i0, i1, i2)
											+ hop_list_.InnerCellIndex(iorbs));
							const integer col =
									memSep
											* (CellIndex(j0, j1, j2)
													+ hop_list_.HopOfOrbs(iorbs)[h].final_idx);
							const scalar val = hop_list_.HopOfOrbs(iorbs)[h].val;
							y[row] += x[col] * val * alpha;
						}
		Hirr_.Multiply(memSep, TotalOfOrbitals(), 1., x, 1., y);
	}

	void ApplyVelocity(const integer memSep, integer dir, scalar* x,
			scalar* y) {
//#pragma omp parallel  for
		for (integer i = 0; i < memSep * TotalOfOrbitals(); i += memSep)
			y[i] = real(0.0) * y[i];

//#pragma omp parallel for 	//This is only for small numbers of cores
		for (int i0 = 0; i0 < CellsInDir(0); i0++)
			for (int i1 = 0; i1 < CellsInDir(1); i1++)
				for (int i2 = 0; i2 < CellsInDir(2); i2++)
					for (int iorbs = 0; iorbs < hop_list_.TotNumRegHops();
							iorbs++)
						for (integer h = 0;
								h < hop_list_.HopOfOrbs(iorbs).size(); h++) {
							const integer j0 = i0
									+ hop_list_.HopOfOrbs(iorbs)[h].Dr[0];
							const integer j1 = i1
									+ hop_list_.HopOfOrbs(iorbs)[h].Dr[1];
							const integer j2 = i2
									+ hop_list_.HopOfOrbs(iorbs)[h].Dr[2];

							const integer row = memSep
									* (CellIndex(i0, i1, i2)
											+ hop_list_.InnerCellIndex(iorbs));
							const integer col =
									memSep
											* (CellIndex(j0, j1, j2)
													+ hop_list_.HopOfOrbs(iorbs)[h].final_idx);
							const scalar val = hop_list_.HopOfOrbs(iorbs)[h].val;
							const real dr = LatticeIndexDifference(i0, i1, i2,
									hop_list_.InnerCellIndex(iorbs), j0, j1, j2,
									hop_list_.HopOfOrbs(iorbs)[h].final_idx,
									dir);
							y[row] += x[col] * val * dr;
						}
	}

	void ApplySpinZ_Velocity(const integer memSep, integer dir, scalar* x,
			scalar* y) {
//#pragma omp parallel  for
		for (integer i = 0; i < memSep * TotalOfOrbitals(); i += memSep)
			y[i] = real(0.0) * y[i];

//#pragma omp parallel for 	//This is only for small numbers of cores
		for (int i0 = 0; i0 < CellsInDir(0); i0++)
			for (int i1 = 0; i1 < CellsInDir(1); i1++)
				for (int i2 = 0; i2 < CellsInDir(2); i2++)
					for (int iorbs = 0; iorbs < hop_list_.TotNumRegHops();
							iorbs++)
						for (integer h = 0;
								h < hop_list_.HopOfOrbs(iorbs).size(); h++) {
							int s0 = hop_list_.InnerCellIndex(iorbs)
									% SpinNumber();
							int sf = hop_list_.HopOfOrbs(iorbs)[h].final_idx
									% SpinNumber();
							if (s0 == sf) {
								const integer j0 = i0
										+ hop_list_.HopOfOrbs(iorbs)[h].Dr[0];
								const integer j1 = i1
										+ hop_list_.HopOfOrbs(iorbs)[h].Dr[1];
								const integer j2 = i2
										+ hop_list_.HopOfOrbs(iorbs)[h].Dr[2];

								const integer row = memSep
										* (CellIndex(i0, i1, i2)
												+ hop_list_.InnerCellIndex(
														iorbs));
								const integer col =
										memSep
												* (CellIndex(j0, j1, j2)
														+ hop_list_.HopOfOrbs(
																iorbs)[h].final_idx);
								const real spin_sgn= 2*s0-1;
								const scalar val =
										hop_list_.HopOfOrbs(iorbs)[h].val*spin_sgn;
								const real dr = LatticeIndexDifference(i0, i1,
										i2, hop_list_.InnerCellIndex(iorbs), j0,
										j1, j2,
										hop_list_.HopOfOrbs(iorbs)[h].final_idx,
										dir);
								y[row] += x[col] * val * dr;
							}
						}
	}

	void ProjKCircle(const std::vector<real> K, const real RK, scalar* x) {

		DFTI_DESCRIPTOR_HANDLE dft_handle;
		///Allocate the memory for the 3D-arra
		scalar *sys_state = new scalar[CellsInDir(2) * CellsInDir(1)
				* CellsInDir(0)];

		//----------------------BEGIN-VALLEY PROJECTION OPERATOR--------------/
		///Pass the information from the linear vector to the 3D array
		for (integer io = 0; io < OrbitalNumber(); ++io)
			for (integer is = 0; is < SpinNumber(); ++is) {
				for (integer i0 = 0; i0 < CellsInDir(0); ++i0)
					for (integer i1 = 0; i1 < CellsInDir(1); ++i1)
						for (integer i2 = 0; i2 < CellsInDir(2); ++i2) {
							const integer k = (CellIndex(i0, i1, i2)
									+ InnerCellIndex(io, is));
							sys_state[i2 * CellsInDir(1) * CellsInDir(0)
									+ i1 * CellsInDir(0) + i0] = x[k];
						}

				mkl::direct_FFT(dft_handle, CellsInDir(2), CellsInDir(1),
						CellsInDir(0), sys_state);

				for (int ik0 = 0; ik0 < CellsInDir(0); ik0++)
					for (int ik1 = 0; ik1 < CellsInDir(1); ik1++) {
						const scalar val = sys_state[ik1 * CellsInDir(0) + ik0];

						//We calculate the k vector
						const real k_state[2] = { ik0 * rec_lat_[0][0]
								/ CellsInDir(0)
								+ ik1 * rec_lat_[1][0] / CellsInDir(1), ik0
								* rec_lat_[0][1] / CellsInDir(0)
								+ ik1 * rec_lat_[1][1] / CellsInDir(1) };

						const real k_dist = sqrt(
								pow(k_state[0] - K[0], 2.0)
										+ pow(k_state[1] - K[1], 2.0));

						//Check if the distance is lesser  and proceed with
						//the calculation
						if (k_dist > RK)
							sys_state[ik1 * CellsInDir(0) + ik0] = 0;
						if (k_dist <= RK) {
							sys_state[ik1 * CellsInDir(0) + ik0] = val;
							//						  NumCal::cout<<k_state[0]<<" "<<k_state[1]<<" "<<RK<<" "<<k_dist <<NumCal::endl;
						}
					}
				///Perform the inverse fourier transform
				mkl::inverse_FFT(dft_handle, CellsInDir(2), CellsInDir(1),
						CellsInDir(0), sys_state);

				///Pass the information back to the linear array
				for (integer i0 = 0; i0 < CellsInDir(0); ++i0)
					for (integer i1 = 0; i1 < CellsInDir(1); ++i1)
						for (integer i2 = 0; i2 < CellsInDir(2); ++i2) {
							const real size = CellsInDir(2) * CellsInDir(1)
									* CellsInDir(0);
							const integer k = (CellIndex(i0, i1, i2)
									+ InnerCellIndex(io, is));
							x[k] = sys_state[i2 * CellsInDir(1) * CellsInDir(0)
									+ i1 * CellsInDir(0) + i0] / size;
						}

			}
		//----------------------END-VALLEY PROJECTION OPERATOR--------------/
		delete[] sys_state;
		DftiFreeDescriptor(&dft_handle);
	}
	;

	void ApplyProjValley_Velocity(const integer memSep, const integer _vidx,
			const real scalR, const integer dir, scalar* x, scalar* y) {

		std::vector<real> K(3, 0);
		switch (_vidx) {
		case -1: {
			K[0] = 2. * ReciprocalLatticeVector(0)[0] / 3.
					+ 1. * ReciprocalLatticeVector(1)[0] / 3.;
			K[1] = 2. * ReciprocalLatticeVector(0)[1] / 3.
					+ 1. * ReciprocalLatticeVector(1)[1] / 3.;
			break;
		}
		case 1: {
			K[0] = 1. * ReciprocalLatticeVector(0)[0] / 3.
					+ 2. * ReciprocalLatticeVector(1)[0] / 3.;
			K[1] = 1. * ReciprocalLatticeVector(0)[1] / 3.
					+ 2. * ReciprocalLatticeVector(1)[1] / 3.;
		}
		}
		//Creates the vector used for the Fourier transform
		scalar* fourier_vec = new scalar[TotalOfOrbitals()];

		const integer Vidx = (_vidx + 1) * 0.5;
		for (integer i0 = 0; i0 < CellsInDir(0); ++i0)
			for (integer i1 = 0; i1 < CellsInDir(1); ++i1)
				for (integer i2 = 0; i2 < CellsInDir(2); ++i2)
					for (integer io = 0; io < OrbitalNumber(); ++io)
						for (integer is = 0; is < SpinNumber(); ++is) {
							const integer k = CellIndex(i0, i1, i2)
									+ InnerCellIndex(io, is);
							fourier_vec[k] = x[memSep * k];
						}
		//This defines the internal reciprocal lattice vectors
		const double Mp[] = { 1. * ReciprocalLatticeVector(0)[0] / 2.
				+ 1. * ReciprocalLatticeVector(1)[0] / 2., 1.
				* ReciprocalLatticeVector(0)[1] / 2.
				+ 1. * ReciprocalLatticeVector(1)[1] / 2. };
		const double RK = scalR
				* sqrt(pow(Mp[0] - K[0], 2.) + pow(Mp[1] - K[1], 2.));

		ProjKCircle(K, RK, fourier_vec);
		for (integer i0 = 0; i0 < CellsInDir(0); ++i0)
			for (integer i1 = 0; i1 < CellsInDir(1); ++i1)
				for (integer i2 = 0; i2 < CellsInDir(2); ++i2)
					for (integer io = 0; io < OrbitalNumber(); ++io)
						for (integer is = 0; is < SpinNumber(); ++is) {
							const integer k = CellIndex(i0, i1, i2)
									+ InnerCellIndex(io, is);
							x[memSep * k] = fourier_vec[k];
						}
		///Apply the velocity  operator to this new x vector
		ApplyVelocity(memSep, dir, x, y);

		for (integer i0 = 0; i0 < CellsInDir(0); ++i0)
			for (integer i1 = 0; i1 < CellsInDir(1); ++i1)
				for (integer i2 = 0; i2 < CellsInDir(2); ++i2)
					for (integer io = 0; io < OrbitalNumber(); ++io)
						for (integer is = 0; is < SpinNumber(); ++is) {
							const integer k = CellIndex(i0, i1, i2)
									+ InnerCellIndex(io, is);
							fourier_vec[k] = y[memSep * k];
						}

		ProjKCircle(K, RK, fourier_vec);

		for (integer i0 = 0; i0 < CellsInDir(0); ++i0)
			for (integer i1 = 0; i1 < CellsInDir(1); ++i1)
				for (integer i2 = 0; i2 < CellsInDir(2); ++i2)
					for (integer io = 0; io < OrbitalNumber(); ++io)
						for (integer is = 0; is < SpinNumber(); ++is) {
							const integer k = CellIndex(i0, i1, i2)
									+ InnerCellIndex(io, is);
							y[memSep * k] = fourier_vec[k];
						}

		delete[] fourier_vec;
	}

	real LatticeIndexDifference(const integer i0, const integer i1,
			const integer i2, const integer iidx, const integer j0,
			const integer j1, const integer j2, const integer jidx,
			const integer iv) {
		real drt = (iidx / SpinNumber() - jidx / SpinNumber()) * Delta[iv];
		real dr[3];

		dr[0] = (i0 - j0) * A[0][iv];
		dr[1] = (i1 - j1) * A[1][iv];
		dr[2] = (i2 - j2) * A[2][iv];
		drt = drt + dr[0] + dr[1] + dr[2];
		return drt;
	}
	;
	void PrintHoppingList() {

		for (int i0 = 0; i0 < 2; i0++)
			for (int i1 = 0; i1 < 2; i1++)
				for (int i2 = 0; i2 < 1; i2++)
					for (int iorbs = 0; iorbs < hop_list_.TotNumRegHops();
							++iorbs)
						for (integer h = 0;
								h < hop_list_.HopOfOrbs(iorbs).size(); h++) {
							NumCal::cout << "H( "
									<< CellIndex(i0, i1, i2)
											+ hop_list_.InnerCellIndex(iorbs)
									<< " , "
									<< CellIndex(
											i0
													+ hop_list_.HopOfOrbs(iorbs)[h].Dr[0],
											i1
													+ hop_list_.HopOfOrbs(iorbs)[h].Dr[1],
											i2
													+ hop_list_.HopOfOrbs(iorbs)[h].Dr[2])
											+ hop_list_.HopOfOrbs(iorbs)[h].final_idx
									<< " ) = "
									<< hop_list_.HopOfOrbs(iorbs)[h].val;
							NumCal::cout << NumCal::endl;
						}

		/*
		 for(  integer h=0;h< hop_array[i].size() ; h++)
		 NumCal::cout<<hop_array[i][h].final_orb<<" "
		 <<hop_array[i][h].final_spin<<" "
		 <<hop_array[i][h].Dr[0]<<" "
		 <<hop_array[i][h].Dr[1]<<" "
		 <<hop_array[i][h].Dr[2]<<" "
		 <<hop_array[i][h].val<<NumCal::endl;
		 */

	}
public:
	std::vector<std::vector<real> > lat_;
	std::vector<std::vector<real> > rec_lat_;
	real volume_;
private:
	std::string label_;
	integer spin_, norb_, coord_;
	std::vector<integer> nsite_;
	bool hamiltonian_setted_;
	IrregularHamiltonian Hirr_;
	//std::vector< real>		lat_;

	HoppingList hop_list_;

};

/** @}*/
}
#endif

