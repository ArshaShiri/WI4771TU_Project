#include <iostream>
#include <typeinfo>
#include <initializer_list>
#include <memory>
#include <cmath>
#include <iomanip>
using namespace std;

class Vector;                                       // is it possible to forward declare the ostream ????
class Matrix;
//////////////////////////////////////////////////////VECTOR CLASS////////////////////////////////////////////////////////////////////////

class Vector
{
public:
	//Constructors
	Vector():length(0),data(0){}
	Vector(int _length):length(_length)
	{
		if (this->length > 0)
		{
            this->data = new double[this->length];
            for (int i=0;i<this->length;i++)                               // If I don't put the date=0, then I will have 
            	this->data[i]=0;                                           // probelms in the mulv functions 
		}

        else
            this->data = 0;
	}

	Vector(initializer_list<double> l):Vector((int)l.size())
    {
    	uninitialized_copy(l.begin(),l.end(),data);
	}

	//Copy constructor
	Vector(const Vector &vector):Vector(vector.length)    //without the copy constructor memory overright in copy operator
    {
        for (int i = 0; i < this->length; i++)             
            data[i] = vector(i);                          // without const double& operator() I'll have an error
    }

    //Destructor
    ~Vector()
    {
    	delete[] data;
    	length=0;
    	//cout<<"Dest  "<<this<<endl;
    }

    //Operation() -> gives the i(th) component of the data.
    double& operator()(int i) { return this->data[i-1]; }
    const double& operator()(int i) const { return this->data[i-1];}  

    //Copy assignment operator 
    Vector& operator=(const Vector &r)
    {
    	if (this != &r)
    	{
    		delete[] data;
    		data=new double [r.length];
    		length=r.length;
    		for (int i=0;i<length;i++)
    			data[i]=r.data[i];
    	}
        return *this;
    }

    //Operation +
    Vector operator+(const Vector &r) const                
    {
        Vector result(this->length);
        for (int i = 0; i < this->length; i++)
            result(i) = this->data[i]+r(i);
        return result;
    }

    //Operation -
    Vector operator-(const Vector &r) const                
    {
        Vector result(this->length);
        for (int i = 0; i < this->length; i++)
            result(i) = this->data[i]-r(i);
        return result;
    }

    //Operation +=
    Vector& operator+=(const Vector &r)                
    {
        for (int i = 0; i < this->length; i++)
            this->data[i]+=r(i);
        return *this;
    }

    //Operation -=
    Vector& operator-=(const Vector &r)                
    {
        for (int i = 0; i < this->length; i++)
            this->data[i]-=r(i);
        return *this;
    }

	//Inner product of 2 vector
    double dot (const Vector &r)
    {
    	double result=0;
    	for (int i = 0; i < this->length; i++)
    		result+=this->data[i]*r(i);
    	return result;
    }

	int length;
private:
	double *data;                                              //why does not work if data is not a pointer 
};

//Scaler vector product
Vector operator*(double l, const Vector &r)                      // why I cannon use const here ????
{
    Vector result(r.length);
    for (int i = 0; i < r.length; i++)
        result(i) = l*r(i);
    return result;
}

//Scaler vector product
Vector operator*(const Vector &r, double l)    
{
    Vector result(r.length);
    for (int i = 0; i < r.length; i++)
        result(i) = l*r(i);
    return result;
}

/*
// ostream overloading for Vector class 
ostream &operator<<(ostream &os, const Vector &v)
{

    for (int i = 0; i < v.length; i++)
    {
    	count = 0
		num = abs(num)
		num = num - int(num)
		while num != 0:
    		num = num * 10
    		count = count + 1
    		num = num - int(num)
    }
    for (int i = 0; i < v.length; i++)
    {
        if (i > 0)
            os << "\n ";
        if (i==0)
        	os <<" ┏"<<setprecision(5)<< v(i)<<" ┓";
        if (i!=0)
        	os <<"┃"<<setprecision(5)<<c v(i)<<"┃";

    }

    return os;
}
*/

ostream &operator<<(ostream &os, const Vector &v)
{
    os << "[";
    for (int i = 0; i < v.length; i++)
    {
        if (i > 0)
            os << ", ";
        os << v(i);
    }
    os << "]";
    return os;
}

//////////////////////////////////////////////////////MATRIX CLASS////////////////////////////////////////////////////////////////////////
class Matrix
{

public:

	//Constructors 
	Matrix():nrows(0),ncols(0),data(0){}
	Matrix(int _n)
	{
		int n=pow(_n,0.5);
		this->nrows=n;
		this->ncols=n;
        this->data = new double[this->nrows*this->ncols];
        for (int i=0;i<this->nrows*this->ncols;i++)
        	this->data[i]=0;

	}

	Matrix(int _nrows,int _ncols):nrows(_nrows),ncols(_ncols)
	{
            this->data = new double[this->nrows*this->ncols];
            for (int i=0;i<this->nrows*this->ncols;i++)
            	this->data[i]=0;
	}

	Matrix(initializer_list<double> l):Matrix((int)l.size())
    {
    	uninitialized_copy(l.begin(),l.end(),data);
	}

	Matrix(int _nrows,int _ncols, initializer_list<double> l):Matrix(_nrows,_ncols)
    {
    	uninitialized_copy(l.begin(),l.end(),data);
	}

	//Copy constructor
	Matrix(const Matrix &matrix):Matrix(matrix.nrows,matrix.ncols)
    {
        for (int i = 0; i < this->nrows*this->ncols; i++)
            (*this)(i) = matrix(i);
    }

    //Destrcutor
	~Matrix()
    {
    	delete[] data;
    	nrows=0;
    	ncols=0;
    }

    //Operation() -> gives the i(th) component of the data.
    double& operator()(int i) { return this->data[i]; }
    const double& operator()(int i) const { return this->data[i];} 
    double& operator()(int i,int j) { return this->data[this->ncols*(i-1)+(j-1)];}
    const double& operator()(int i,int j) const { return this->data[this->ncols*(i-1)+(j-1)];} 

    //Copy operator
    Matrix& operator=(const Matrix &r)
    {
    	if (this != &r)
    	{
    		delete[] data;
    		data=new double [r.nrows*r.ncols];
    		nrows=r.nrows;
    		ncols=r.ncols;
    		for (int i=0;i<nrows*ncols;i++)
    			data[i]=r.data[i];
    	}
        return *this;
    }

    //Operation +
    Matrix operator+(const Vector &r) const                // I get a warning when I use Matrix& why ??
    {
        Matrix result(this->nrows,this->ncols);
        for (int i = 0; i < this->nrows*this->ncols; i++)
            result(i) = this->data[i]+r(i);
        return result;
    }

    //Operation -
    Matrix operator-(const Vector &r) const                
    {
        Matrix result(this->nrows,this->ncols);
        for (int i = 0; i < this->nrows*this->ncols; i++)
            result(i) = this->data[i]-r(i);
        return result;
    }

    //Operation +=
    Matrix& operator+=(const Vector &r)             
    {
        for (int i = 0; i < this->nrows*this->ncols; i++)
            this->data[i]+=r(i);
        return *this;
    }

    //Operation -=
    Matrix& operator-=(const Vector &r)              
    {
        for (int i = 0; i < this->nrows*this->ncols; i++)
            this->data[i]-=r(i);
        return *this;
    }

	//template<typename A,typename B> 
    Matrix mul(const Matrix &m) const
    {
    	Matrix result(nrows,m.ncols);
    	for (int k=0;k < nrows; k++)
    	{
    		for (int i = 0; i < m.ncols; i++)
    		{
    			for (int j = 0; j < m.nrows; j++)
    				result(k*m.ncols+i)+=this->data[k*ncols+j]*m(j*m.ncols+i);
    		}
    	}	
    	return result;
    }

    Vector vmul(const Vector &v) const
    {
    	Vector result(ncols);
    	for (int k=0;k < ncols; k++)
    	{
    		for (int i = 0; i < v.length; i++)
    			result(k)+=this->data[i*this->ncols+k]*v(i);	
    	}	
    	return result;
    }

    Vector mulv(const Vector &v) const
    {
    	Vector result(nrows);
    	for (int k=0;k < nrows; k++)
    	{
    		for (int i = 0; i < nrows; i++)
    		{
    			result(k)+=this->data[k*ncols+i]*v(i);
    		}

    	}	
    	return result;	
    }

    Matrix transpose()
    {
    	Matrix dummy(this->ncols,this->nrows);
    	for (int i=0;i<dummy.nrows;i++)
    	{
    			for (int j=0;j<dummy.ncols;j++)
  					dummy(i*dummy.ncols+j)=this->data[j*this->ncols+i];
    	}
    	return dummy;
    }

    Vector CG(const Vector &v)
    {
    	Vector result(v.length);
    	for (int i=0;i<v.length;i++)
    		result(i)=i;
    	double alfa;
    	double beta;
    	Vector r(v.length);
    	Vector p(v.length);
    	r=v-this->mulv(result);
    	p=r;
    	double rold=r.dot(r);
    	for (int j=0;j<v.length;j++)
    	{   		
    		
    		alfa=(rold)/(this->vmul(p).dot(p));
    		result+=alfa*p;
    		r=r-alfa*(this->mulv(p));
    		double rnew=(r.dot(r));
    		p=r+(rnew/rold)*p;
    		rold=rnew;	
    	}
    	return result;
    }

	int nrows;
	int ncols;
private:
	double* data;
};


Matrix operator*(double l, const Matrix &r)    
{
    Matrix result(r.nrows,r.ncols);
    for (int i = 0; i < r.nrows*r.ncols; i++)
        result(i) = l*r(i);
    return result;
}

Matrix operator*(const Matrix &r, double l)    
{
    Matrix result(r.nrows,r.ncols);
    for (int i = 0; i < r.nrows*r.ncols; i++)
        result(i) = l*r(i);
    return result;
}

    //Scaler matrix product
	//Matrix operator*(double l)
	//{
    //	Matrix result(nrows,ncols);
    //	for (int i = 0; i < nrows*ncols; i++)
      //  	result(i) = l*this->data[i];
    	//return result;
	//}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ostream &operator<<(ostream &os, const Matrix &m)
{
	os<<endl;
	os << "{";
	for (int j=0; j<m.nrows;j++)
	{
		if (j > 0)
		 os << ", \n";
    	os << "[";
    	for (int i = j*m.ncols; i <j*m.ncols+m.ncols; i++)
    	{
        	if (i != j*m.ncols)
            	os << ", ";
        	os << m(i);
   		}
    	os << "]";
    }
    os << "}";
    return os;
}