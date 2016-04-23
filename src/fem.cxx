#include "integration.hpp"
#include "vector.hpp"
#include <iostream>
#include <functional>

class NodalCoordinates: public Matrix
{
public:
	using Matrix::Matrix;
};

class Connect:public Matrix
{
public:
	using Matrix::Matrix;
};

class BC:public Matrix
{
public:
	using Matrix::Matrix;
};

class Tri3
{
public:		//why can't I put Matrix K(6,6) as an atrribute ??
	Tri3(double e,double nu):E(e),Nu(nu){}
 														//check this with matlab code

	Matrix D()
	{	
		double lambda=(Nu*E)/((1+Nu)*(1-2*Nu));
		double mu=E/(2*(1+Nu));
		Matrix d(3,3,{lambda+2*mu,lambda,0.0,lambda,lambda+2*mu,0.0,0.0,0.0,mu});
		return d;
	}

	Matrix B(Matrix C)
	{
		Matrix j(2,2);
		Matrix dn(2,3,{-1.0,1.0,0.0,-1.0,0.0,1.0});
		j=dn.mul(C.transpose());
		double detj=(j(1,1)*j(2,2))-(j(1,2)*j(2,1));
		Matrix invj(2,2,{j(2,2),-j(1,2),-j(2,1),j(1,1)});
		invj=1/detj*invj;
		dn=invj.mul(dn);
		Matrix b(3,6,{dn(1,1),0.0,dn(1,2),0.0,dn(1,3),0.0,0.0,dn(2,1),0.0,dn(2,2),0.0,dn(2,3),dn(2,1),dn(1,1),dn(2,2),dn(1,2),dn(2,3),dn(1,3)});
		return b;
	}

	Matrix K(Matrix C)
	{
		Matrix k(6,6);
		k=(((this->B(C)).transpose()).mul(this->D())).mul(this->B(C));
		Gauss<6,double> a(1);
		Dim<6,double> f;
		f.coordinate[0]=C(1,1);
		f.coordinate[1]=C(2,1);
		f.coordinate[2]=C(1,2);
		f.coordinate[3]=C(2,2);
		f.coordinate[4]=C(1,3);
		f.coordinate[5]=C(2,3);
		auto g = [](double x,double y) { return 1; };
		k=k*a.integrate(g,f);
		return k;
	}

	Matrix Solve(NodalCoordinates N,Connect C)
	{
		Matrix GK(N.ncols,N.ncols);
		for (auto i=0;i<C.nrows;i++)
		{
			Matrix LK(6,6);
			Matrix el(2,3);
			for (auto j=0;j<C.ncols;j++)
			{
				el(1,j+1)=N(C(i+1,j+1),1);
				el(2,j+1)=N(C(i+1,j+1),2);
				LK=this->K(el);
				
			} 
			cout<<el<<endl;
			cout<<LK<<endl;
		}
		return GK;
	}

	double E;
	double Nu;
};

int main()
{
	/*
	NodalCoordinates A(2,3,{0.0,1.0,0.0,0.0,0.0,1.0});
	Tri3 B(1,0.2);
	cout<<B.B(A)<<endl;
	cout<<B.D()<<endl;
	cout<<B.K(A)<<endl;

	NodalCoordinates C(2,3,{0.0,3.0,3.0,0.0,0.0,3.0});
	Tri3 D(1,0.27);
	cout<<D.B(C)<<endl;
	cout<<D.K(C)<<endl;

	NodalCoordinates E(2,3,{0.0,3.0,0.0,0.0,3.0,3.0});
	Tri3 F(1,0.27);
	cout<<F.B(E)<<endl;
	cout<<F.K(E)<<endl;
	*/
	NodalCoordinates A(4,2,{0.0,0.0,3.0,0.0,3.0,3.0,0.0,3.0});
	Connect B(2,3,{1.0,2.0,3.0,1.0,3.0,4.0});
	Tri3 C(1,0.27);
	C.Solve(A,B);

    return 0;
}
