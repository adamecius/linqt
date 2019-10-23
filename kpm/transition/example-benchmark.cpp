/*
 * example-benchmark.cpp
 *
 *  Created on: 24/08/2016
 *      Author: jgarcia
 */
#include "types_definitions.hpp"
#include "kpm.hpp"
#include "sys/time.h"

int main()
{


	///Dimension of the system
	const my::integer dim=2*1000*1000;
	///Number of moments
	const my::integer M=10;
	///The Kpm Class
	Kpm	kpm_example("conductividad",dim);
	///Print the KPM usage
	kpm_example.PrintKPMLoad();


	///Measure of time
	struct timeval start, end;


	std::cout<<"\nIn order to check the linear algebra functions"<<std::endl
			 <<"First we create a random vector"<<std::endl;
	kpm_example.RandomFill(0,0);
	if(dim<10)
		kpm_example.PrintVector(0,0);


	std::cout<<"The first test is to compute the norm:"
			 <<" using the standandr nrm2 function:"<<std::endl;
	const my::real
	norm_ver1=kpm_example.kpm_linalg.nrm2(kpm_example.ChebVec(0,0));
	std::cout<<norm_ver1<<std::endl;

	std::cout<<"using the dot product:"<<std::endl;
	const my::real
	norm_ver2=sqrt(kpm_example.kpm_linalg.dot(kpm_example.ChebVec(0,0),kpm_example.ChebVec(0,0)).real());
	std::cout<<norm_ver2<<std::endl;


	std::cout<<"The normalize vector is:"<<std::endl;
	kpm_example.kpm_linalg.normalize(kpm_example.ChebVec(0,0));
	if(dim<10)
		kpm_example.PrintVector(0,0);

	std::cout<<"and effectively the norm is:"<<std::endl;
	const my::real
	norm_ver3=kpm_example.kpm_linalg.nrm2(kpm_example.ChebVec(0,0));
	std::cout<<norm_ver3<<std::endl;


	const my::real alpha=2;
	std::cout<<"If we scale the vector by: "<<alpha<<std::endl;
	std::cout<<"Then the vector becomes"<<std::endl;
	kpm_example.kpm_linalg.scale(alpha,kpm_example.ChebVec(0,0));
	if(dim<10)
		kpm_example.PrintVector(0,0);
	std::cout<<"which has norm:"<<std::endl;
	const my::real
	newnorm=kpm_example.kpm_linalg.nrm2(kpm_example.ChebVec(0,0));
	std::cout<<newnorm<<std::endl;


	std::cout<<"A duplication can be also done by copying the vector"<<alpha<<" :"<<std::endl;
	kpm_example.kpm_linalg.copy(kpm_example.ChebVec(0,0),kpm_example.ChebVec(1,0));
	const my::real a=1;
	std::cout<<"then set a="<<a<<" and use the axpy function"<<std::endl;
	kpm_example.kpm_linalg.axpy(a,kpm_example.ChebVec(0,0),kpm_example.ChebVec(1,0));
	std::cout<<"which lead to the vector: "<<std::endl;
	if(dim<10)
		kpm_example.PrintVector(1,0);
	std::cout<<"which has norm:"<<std::endl;
	const my::real
	newnorm2=kpm_example.kpm_linalg.nrm2(kpm_example.ChebVec(1,0));
	std::cout<<newnorm2<<std::endl;

	/*
	kpm_example.kpm_linalg.copy(kpm_example.ChebVec(0,0),kpm_example.ChebVec(1,0));
	my::scalar norm1,norm2;

	/*

	gettimeofday(&start, NULL);


	for(int m=0;m<M;m++)
		for(int n=0;n<M;n++)
		{
			//		kpm_example.ChebyshevIteration( 0, m );
			//		kpm_example.ChebyshevIteration( 1, n );
			kpm_example.kpm_linalg.copy(kpm_example.ChebVec(0,0),kpm_example.ChebVec(0,1));
			kpm_example.kpm_linalg.copy(kpm_example.ChebVec(0,0),kpm_example.ChebVec(0,2));
			kpm_example.kpm_linalg.copy(kpm_example.ChebVec(1,0),kpm_example.ChebVec(1,1));
			kpm_example.kpm_linalg.copy(kpm_example.ChebVec(1,0),kpm_example.ChebVec(1,2));
		}
	norm1=kpm_example.kpm_linalg.dot(kpm_example.ChebVec(0,0),kpm_example.ChebVec(0,0));
	std::cout<<"Unormalized norm "<<norm1<<std::endl;

	kpm_example.kpm_linalg.normalize(kpm_example.ChebVec(0,0));

	norm1=kpm_example.kpm_linalg.dot(kpm_example.ChebVec(0,0),kpm_example.ChebVec(0,0));
	std::cout<<"normalized norm "<<norm1<<std::endl;
	gettimeofday(&end, NULL);

//	std::cout<<"The random vector was filled in :"<< (end.tv_sec - start.tv_sec)*1000000 + (end.tv_usec - start.tv_usec)<<" ms"<<std::endl;
	//kpm_example.PrintVector(0,0);
*/

	return 0;
}



