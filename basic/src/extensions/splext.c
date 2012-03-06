#include <Python.h>
#include "numpy/oldnumeric.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

static PyObject
*splext_secderivs(PyObject *self, PyObject *args)
{
	PyArrayObject *x;
	PyArrayObject *y;
	PyArrayObject *secder;
	double *y2,*u;
	double tmp;
	double p,qn,un,sig;
	int nvalues;
	int xs,ys;
	char *xd,*yd;
	int dims[1];
	int i;

	//check Arguments
	if(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &x, &PyArray_Type,&y))
		return NULL;
	if (x->descr->type_num != PyArray_DOUBLE || x->nd != 1 ||
		y->descr->type_num != PyArray_DOUBLE || y->nd != 1 || x->dimensions[0]!=y->dimensions[0])
	{
		PyErr_SetString(PyExc_ValueError, "x and y must be 1 dimensional arrays of type float and equal length.");
		return NULL;
	}
	nvalues=x->dimensions[0];
	//store data array bases and strides in temp variables
	//to make code shorter and more readable
	//(should also save some memory accesses)
	xd=x->data;
	yd=y->data;
	xs=x->strides[0];
	ys=y->strides[0];
	//allocate temporary arrays
	y2=malloc(nvalues*sizeof(double));
	u=malloc(nvalues*sizeof(double));

	//set natural spline boundary conditions
	y2[0]=0.0;
	u[0]=0.0;
	y2[nvalues-1]=0.0;
	u[nvalues-1]=0.0;
	//initialize return PyArrayObject
	dims[0]=nvalues;
	secder = (PyArrayObject *)PyArray_FromDims(1,dims,PyArray_DOUBLE);
	//first loop of tridiagonal algorithm
	for (i=1;i<nvalues-1;i++)
	{
		sig=(*(double *)(xd+i*xs)-*(double *)(xd+(i-1)*xs))/(*(double *)(xd+xs*(i+1))- *(double *)(xd+xs*(i-1)));
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		tmp=(*(double *)(yd+ys*(i+1))- *(double *)(yd+ys*i))/(*(double *)(xd+xs*(i+1))- *(double *)(xd+xs*i));
		tmp-=(*(double *)(yd+ys*i)- *(double *)(yd+ys*(i-1)))/(*(double *)(xd+xs*i)- *(double *)(xd+xs*(i-1)));
		u[i]=(6.0*tmp/(*(double *)(xd+xs*(i+1))- *(double *)(xd+xs*(i-1)))-sig*u[i-1])/p;
		}

	//"backsubstitution" loop of tridiagonal algorithm
	for (i=nvalues-2;i>=0;i--)
	{
		y2[i]=y2[i]*y2[i+1]+u[i];
		*(double *)(secder->data + i*secder->strides[0])=y2[i];
	}
	//free temporary arrays and return
	free(u);
	free(y2);
	return PyArray_Return(secder);
}


static PyObject
*splext_splint(PyObject *self, PyObject *args)
{
	PyArrayObject *x;
	PyArrayObject *y;
	PyArrayObject *secder;
	double xint;

	int hi,lo,i,nvalues;
	double dx,a,b,yout,tmp;

	//check arguments
	if(!PyArg_ParseTuple(args, "O!O!O!d", &PyArray_Type, &x, &PyArray_Type,&y, &PyArray_Type, &secder, &xint))
		return NULL;
	if (x->descr->type_num != PyArray_DOUBLE || x->nd != 1 ||
		y->descr->type_num != PyArray_DOUBLE || y->nd != 1 ||
		secder->descr->type_num != PyArray_DOUBLE || secder->nd != 1 
		|| x->dimensions[0]!=y->dimensions[0] ||
		x->dimensions[0]!=secder->dimensions[0])
	{
		PyErr_SetString(PyExc_ValueError, "x,y and secder must be 1 dimensional arrays of type float and equal length.");
		return NULL;
	}
	//this might be redundant, but out-of-fitrange interpolation is generally a nono
	if ((xint < *(double *)(x->data)) || xint > *(double *)(x->data + x->strides[0] * (x->dimensions[0]-1)))
	{
		PyErr_SetString(PyExc_ValueError, "x value out of fit data range");
		return NULL;
	}

	nvalues=x->dimensions[0];

	//find x-interval in which to interpolate
	for (i=1;(*(double *)(x->data + x->strides[0]*i) < xint) && (i < nvalues-1);i++);
	hi=i;
	lo=i-1;
	
	//cspline evaluation
	dx=*(double *)(x->data + x->strides[0]*hi)- *(double *)(x->data + x->strides[0] * lo);
	a=(*(double *)(x->data + x->strides[0]*hi) - xint)/dx;
	b=1-a;
	yout= a* *(double *)(y->data + y->strides[0]*lo);
	yout+=b* *(double *)(y->data + y->strides[0]*hi);
	tmp= ((a*a*a-a)* *(double *)(secder->data + secder->strides[0] * lo));
	tmp+=((b*b*b-b)* *(double *)(secder->data + secder->strides[0] * hi));
	tmp*=(dx*dx)/6.0;
	yout+=tmp;

	//return
	return Py_BuildValue("d", yout);
}



static PyObject
*splext_splder(PyObject *self, PyObject *args)
{
	PyArrayObject *x;
	PyArrayObject *y;
	PyArrayObject *secder;
	double xint;

	int hi,lo,i,nvalues;
	double dx,a,b,derout,tmp;

	//parse and check arguments
	if(!PyArg_ParseTuple(args, "O!O!O!d", &PyArray_Type, &x, &PyArray_Type,&y, &PyArray_Type, &secder, &xint))
		return NULL;
	if (x->descr->type_num != PyArray_DOUBLE || x->nd != 1 ||
		y->descr->type_num != PyArray_DOUBLE || y->nd != 1 ||
		secder->descr->type_num != PyArray_DOUBLE || secder->nd != 1 
		|| x->dimensions[0]!=y->dimensions[0] ||
		x->dimensions[0]!=secder->dimensions[0])
	{
		PyErr_SetString(PyExc_ValueError, "x,y and secder must be 1 dimensional arrays of type float and equal length.");
		return NULL;
	}
	//might be redundant, see splext_splint
	if ((xint < *(double *)(x->data)) || xint > *(double *)(x->data + x->strides[0] * (x->dimensions[0]-1)))
	{
		PyErr_SetString(PyExc_ValueError, "x value out of fit data range");
		return NULL;
	}

	nvalues=x->dimensions[0];
	
	//find x range in which to compute 1st derivative
	for (i=1;(*(double *)(x->data + x->strides[0]*i) < xint) && (i < nvalues-1);i++);
	hi=i-0;
	lo=i-1;
	
	dx=*(double *)(x->data + x->strides[0]*hi)- *(double *)(x->data + x->strides[0] * lo);
	a=(*(double *)(x->data + x->strides[0]*hi) - xint)/dx;
	b=1.0-a;
	derout=(*(double *)(y->data + y->strides[0]*hi)- *(double *)(y->data+y->strides[0]*lo))/dx;
	tmp= -((3.0*a*a)-1.0)* *(double *)(secder->data+secder->strides[0]*lo);
	tmp+=(( 3.0*b*b)-1.0)* *(double *)(secder->data+secder->strides[0]*hi);
	tmp*= (dx)/6.0;
	derout+=tmp;
	
	//return
	return Py_BuildValue("d", derout);
}




static PyMethodDef SplExtMethods[] = {
	{"spline", splext_secderivs, METH_VARARGS,
		"return cubic spline second derivatives array for x and y float value grid arrays"},
		{"splint", splext_splint, METH_VARARGS,
		"return cspline interpolation of the cspline specified by x,y,secder at point xint"},
		{"splder", splext_splder, METH_VARARGS,
		"return cspline derivative of the cspline specified by x,y,secder at point xint"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initsplext(void)
{
    (void) Py_InitModule("splext", SplExtMethods);
	import_array();
}
