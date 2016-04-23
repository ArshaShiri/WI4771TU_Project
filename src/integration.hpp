#include <iostream>
#include <functional>
#include <cmath>
#include <array>

using namespace std;

template<int dim, typename T>
class Dim
{
public:
    // constructors
    Dim()
    {
        for (auto i=0; i<dim; i++)
            coordinate[i] = 0;
    }
    
    Dim(std::array<T,dim> c)
    {
        for (auto i=0; i<dim; i++)
            coordinate[i] = c[i];
    }

    std::array<T,dim> coordinate;
};


template<int dim, typename T>
class Quadrature
{
public:
    Quadrature():numpoints(0),weights(nullptr),intpoints(nullptr){}
    
    Quadrature(int n):numpoints(n),weights(new double[n]),intpoints(new double[n]){}
          
    ~Quadrature()
    {
        delete[] weights;
        delete[] intpoints;
        numpoints = 0;
    }

    virtual double jacobian(int n,int q,Dim<dim,T> a) = 0;
    
    virtual double mapping1d(int n,Dim<dim,T> a) = 0;
    
    virtual Dim<2,T> mapping2d(int a,int b,Dim<dim,T> c) = 0;
    
    virtual T integrate(const double (*func)(const T x), Dim<dim,T> a)
    {
        T integral(0);
        for (auto i=0; i<numpoints; i++)
            integral+=weights[i]*func(mapping1d(i,a));
        return jacobian(1,1,a)*integral;
    }

    // Virtual integrator: lambda expression
    virtual T integrate(function<T(T)> func, Dim<dim,T> a)
    {
        T integral(0);
        for (auto i=0; i<numpoints; i++)
            integral += this->weights[i]*func(mapping1d(i,a));
        return jacobian(1,1,a)*integral;
    }
    
protected: 																	
    double* weights;
    double* intpoints;
    int numpoints;
};


//Gaussian integration in 1d
template<int dim, typename T>
class Gauss : public Quadrature<dim, T>
{
public:
    
    virtual Dim<2,T> mapping2d(int a,int b,Dim<dim,T> c){return Dim<2,double>();} 

    Gauss(int n):Quadrature<dim,T>(n)
    {
        switch(this->numpoints)
        {
        case 1:
            this->weights[0]={2.0};
            this->intpoints[0]={0.0};
            break;
        case 2:
            this->weights[0]={1.0};
            this->weights[1]={1.0};
            this->intpoints[0]={-0.57735026919};
            this->intpoints[1]={0.57735026919};
            break;
        case 3:
            this->weights[0]={5.0/9.0};
            this->weights[1]={8.0/9.0};
            this->weights[2]={5.0/9.0};
            this->intpoints[0]={-0.774596669241};
            this->intpoints[1]={0.0};
            this->intpoints[2]={0.774596669241};
            break;
        default:
            cout << "Not defined" <<endl;
            exit(1);
        }
    }

    virtual double jacobian(int xx,int yy,Dim<dim,T> a)							
    {
        return (a.coordinate[1]-a.coordinate[0])/2.0;
    }
    
    virtual double mapping1d(int n,Dim<dim,T> a)
    {
        double mapping(0);
            mapping=(a.coordinate[0]*((1-this->intpoints[n])/2.0)+a.coordinate[1]*((1+this->intpoints[n])/2.0));
        return mapping;
    }
};


//Gaussian integration in 2d for rectangle.
template<typename T>                                                         //why should I put <8,T> for Guass here ????
class Gauss<8,T> : public Quadrature<8,T> 
{

public:
    
    virtual double mapping1d(int n,Dim<8,T> a){return 1;}
        
    Gauss(int n):Quadrature<8,T>(n)
    {
        switch(this->numpoints)
        {
        case 1:
            this->weights[0]={2.0};
            this->intpoints[0]={0.0};
            break;
        case 2:
            this->weights[0]={1.0};
            this->weights[1]={1.0};
            this->intpoints[0]={-0.57735026919};
            this->intpoints[1]={0.57735026919};
            break;
        case 3:
            this->weights[0]={5.0/9.0};
            this->weights[1]={8.0/9.0};
            this->weights[2]={5.0/9.0};
            this->intpoints[0]={-0.774596669241};
            this->intpoints[1]={0.0};
            this->intpoints[2]={0.774596669241};
            break;
        default:
            cout << "Not defined" <<endl;
            exit(1);
        }
    } 

    virtual double jacobian(int xx,int yy,Dim<8,T> a)
    {   
        double dxN1= -1.0/4.0*(1-this->intpoints[yy]);
        double dxN2= 1.0/4.0*(1-this->intpoints[yy]);
        double dxN3= 1.0/4.0*(1+this->intpoints[yy]);
        double dxN4= -1.0/4.0*(1+this->intpoints[yy]);

        double dyN1= -1.0/4.0*(1-this->intpoints[xx]);
        double dyN2= -1.0/4.0*(1+this->intpoints[xx]);
        double dyN3= 1.0/4.0*(1+this->intpoints[xx]);
        double dyN4= 1.0/4.0*(1-this->intpoints[xx]);

        //double result(dxN1); 
        double result=(dxN1*a.coordinate[0]+dxN2*a.coordinate[2]+dxN3*a.coordinate[4]+dxN4*a.coordinate[6])*
        (dyN1*a.coordinate[1]+dyN2*a.coordinate[3]+dyN3*a.coordinate[5]+dyN4*a.coordinate[7])-
        (dyN1*a.coordinate[0]+dyN2*a.coordinate[2]+dyN3*a.coordinate[4]+dyN4*a.coordinate[6])*
        (dxN1*a.coordinate[1]+dxN2*a.coordinate[3]+dxN3*a.coordinate[5]+dxN4*a.coordinate[7]);
        return result;
    }


    virtual Dim<2,T> mapping2d(int xx,int yy,Dim<8,T> a)
    {   
        Dim<2,T> result;
        double N1=1.0/4.0*(1-this->intpoints[xx])*(1-this->intpoints[yy]);
        double N2=1.0/4.0*(1+this->intpoints[xx])*(1-this->intpoints[yy]);
        double N3=1.0/4.0*(1+this->intpoints[xx])*(1+this->intpoints[yy]);
        double N4=1.0/4.0*(1-this->intpoints[xx])*(1+this->intpoints[yy]);
        result.coordinate[0]=N1*a.coordinate[0]+N2*a.coordinate[2]+N3*a.coordinate[4]+N4*a.coordinate[6];
        result.coordinate[1]=N1*a.coordinate[1]+N2*a.coordinate[3]+N3*a.coordinate[5]+N4*a.coordinate[7];
        return result;
    }

    virtual T integrate(const T (*func)(const T x,T y),Dim<8,T> a)
    {
        T integral(0); 
        for (auto i=0; i<this->numpoints; i++)
        {
            for (auto j=0; j<this->numpoints; j++)
            integral+=jacobian(i,j,a)*
        this->weights[j]*this->weights[i]*func(mapping2d(i,j,a).coordinate[0],mapping2d(i,j,a).coordinate[1]);
        }
        return integral;
    } 

    virtual T integrate(function<T(T,T)> func, Dim<8,T> a)
    {
        T integral(0); 
        for (auto i=0; i<this->numpoints; i++)
        {
            for (auto j=0; j<this->numpoints; j++)
            integral+=jacobian(i,j,a)*
        this->weights[j]*this->weights[i]*func(mapping2d(i,j,a).coordinate[0],mapping2d(i,j,a).coordinate[1]);
        }
        return integral;
    }

};



template<typename T>                                                         
class Gauss<6,T> : public Quadrature<6,T>
{
public:	

	virtual double mapping1d(int n,Dim<6,T> a){return 1;}
	
    Gauss(int n):Quadrature<6,T>(n)
    {
        switch(this->numpoints)
        {
        case 1:
            this->weights[0]={1.0/2.0};
            this->intpoints[0]={1.0/3.0};
            break;
        case 3:
            this->weights[0]={1.0/6.0};
            this->weights[1]={1.0/6.0};
            this->weights[2]={1.0/6.0};
            this->intpoints[0]={0};
            this->intpoints[1]={1.0/2.0};
            this->intpoints[2]={1.0/2.0};
            break;
        default:
            cout << "Not defined" <<endl;
            exit(1);
        }
    } 

    //Guass2D():Gauss1D(n){}
    virtual double jacobian(int xx,int yy,Dim<6,T> a)
    {	
    	double dxN1= -1.0;
    	double dxN2= 1.0;
    	double dxN3= 0.0;

    	double dyN1= -1.0;
    	double dyN2= 0.0;
    	double dyN3= 1.0;
    	
    	 //double result=
    	double result=(dxN1*a.coordinate[0]+dxN2*a.coordinate[2]+dxN3*a.coordinate[4])*
    	(dyN1*a.coordinate[1]+dyN2*a.coordinate[3]+dyN3*a.coordinate[5])-
    	(dyN1*a.coordinate[0]+dyN2*a.coordinate[2]+dyN3*a.coordinate[4])*
    	(dxN1*a.coordinate[1]+dxN2*a.coordinate[3]+dxN3*a.coordinate[5]);
    	return result;
    }
    
    virtual Dim<2,T> mapping2d(int xx,int yy,Dim<6,T> a)
    {
    	Dim<2,T> result;
    	double N1=1-this->intpoints[xx]-this->intpoints[yy];
    	double N2=this->intpoints[xx];
    	double N3=this->intpoints[yy];
    	result.coordinate[0]=N1*a.coordinate[0]+N2*a.coordinate[2]+N3*a.coordinate[4];
    	result.coordinate[1]=N1*a.coordinate[1]+N2*a.coordinate[3]+N3*a.coordinate[5];
    	return result;
    }

    virtual T integrate(const T (*func)(const T x,T y),Dim<6,T> a)
    {
        T integral(0);
        Dim<3,int> j({1,0,1}); 
        for (auto i=0; i<this->numpoints; i++)
        {
            integral+=jacobian(i,j.coordinate[i],a)*
        this->weights[i]*func(mapping2d(i,j.coordinate[i],a).coordinate[0],mapping2d(i,j.coordinate[i],a).coordinate[1]);
        }
        return integral;
    }

    virtual T integrate(function<T(T,T)> func, Dim<6,T> a)
    {
        T integral(0);
        Dim<3,int> j({1,0,1}); 
        for (auto i=0; i<this->numpoints; i++)
        {
            integral+=jacobian(i,j.coordinate[i],a)*
        this->weights[i]*func(mapping2d(i,j.coordinate[i],a).coordinate[0],mapping2d(i,j.coordinate[i],a).coordinate[1]);
        }
        return integral;
    }
};