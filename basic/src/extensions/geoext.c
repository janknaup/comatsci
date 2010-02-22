#include <Python.h>
#include "numpy/oldnumeric.h"

#include <math.h>

static PyObject *
geoext_dmatrix(PyObject *self, PyObject *args)
{
	PyArrayObject *pos;
	PyArrayObject *dmatrix;
	int natoms;
	int i,j;
	double dx,dy,dz,distance;
	int dimensions[2];
	
	
	if(!PyArg_ParseTuple(args, "O!", &PyArray_Type, &pos))
		return NULL;
	if (pos->nd !=2 || pos->descr->type_num != PyArray_DOUBLE
			|| pos->dimensions[1]!=3){
		PyErr_SetString(PyExc_ValueError, "array must be two-dimensional, of shape (n,3) and of type float");
		return NULL;
	}
	natoms=pos->dimensions[0];
	dimensions[0]=natoms;
	dimensions[1]=natoms;
	dmatrix = (PyArrayObject *)PyArray_FromDims(2,dimensions,PyArray_DOUBLE);
	for (i=0;i<natoms;i++)
	{
		for (j=i;j<natoms;j++)
		{
			if (j!=i)
			{
				dx=*(double *)(pos->data + i * pos->strides[0])
					- *(double *)(pos->data + j * pos->strides[0]);
				dy=*(double *)(pos->data + i * pos->strides[0] + pos->strides[1])
					- *(double *)(pos->data + j * pos->strides[0] + pos->strides[1]);
				dz=*(double *)(pos->data + i * pos->strides[0] + 2 * pos->strides[1])
					- *(double *)(pos->data + j * pos->strides[0] + 2 * pos->strides[1]);
				distance=sqrt((dx*dx)+(dy*dy)+(dz*dz));
				*(double *)(dmatrix->data + i*dmatrix->strides[0]+j*dmatrix->strides[1])=distance;
				*(double *)(dmatrix->data + j*dmatrix->strides[0]+i*dmatrix->strides[1])=distance;
			} else {
				*(double *)(dmatrix->data + i*dmatrix->strides[0]+j*dmatrix->strides[1])=0.0;
			}
		}
	}
	return PyArray_Return(dmatrix);
}
	

static PyObject
*geoext_sdmatrix(PyObject *self, PyObject *args)
{
	PyArrayObject *pos;
	PyArrayObject *dmatrix;
	PyArrayObject *lattice;
	int natoms;
	int i,j,k,l,m;
	double dx,dy,dz,distance;
	double pdx,pdy,pdz;
	double ax,ay,az,bx,by,bz,cx,cy,cz;
	int dimensions[2];
	
	
	if(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &pos, &PyArray_Type,&lattice))
		return NULL;
	if (pos->nd !=2 || pos->descr->type_num != PyArray_DOUBLE
			|| pos->dimensions[1]!=3){
		PyErr_SetString(PyExc_ValueError, "geometry array must be two-dimensional, of shape (n,3) and of type float");
		return NULL;
	}
	if (lattice->nd !=2 || lattice->descr->type_num != PyArray_DOUBLE
			|| lattice->dimensions[0]!=3 || lattice->dimensions[1] !=3){
		PyErr_SetString(PyExc_ValueError, "lattice array must be two-dimensional, of shape (3,3) and of type float");
		return NULL;
	}
	
	natoms=pos->dimensions[0];
	dimensions[0]=natoms;
	dimensions[1]=natoms;
	dmatrix = (PyArrayObject *)PyArray_FromDims(2,dimensions,PyArray_DOUBLE);
	ax= *(double *)(lattice->data);
	ay= *(double *)(lattice->data+lattice->strides[1]);
	az= *(double *)(lattice->data+2*lattice->strides[1]);
	bx= *(double *)(lattice->data+lattice->strides[0]);
	by= *(double *)(lattice->data+lattice->strides[1]+lattice->strides[0]);
	bz= *(double *)(lattice->data+2*lattice->strides[1]+lattice->strides[0]);
	cx= *(double *)(lattice->data+2*lattice->strides[0]);
	cy= *(double *)(lattice->data+lattice->strides[1]+2*lattice->strides[0]);
	cz= *(double *)(lattice->data+2*lattice->strides[1]+2*lattice->strides[0]);
	for (i=0;i<natoms;i++)
	{
		for (j=i;j<natoms;j++)
		{
			dx=*(double *)(pos->data + i * pos->strides[0])
				- *(double *)(pos->data + j * pos->strides[0]);
			dy=*(double *)(pos->data + i * pos->strides[0] + pos->strides[1])
				- *(double *)(pos->data + j * pos->strides[0] + pos->strides[1]);
			dz=*(double *)(pos->data + i * pos->strides[0] + 2 * pos->strides[1])
				- *(double *)(pos->data + j * pos->strides[0] + 2 * pos->strides[1]);
			*(double *)(dmatrix->data + i*dmatrix->strides[0]+j*dmatrix->strides[1])=1000000.0;
			*(double *)(dmatrix->data + j*dmatrix->strides[0]+i*dmatrix->strides[1])=1000000.0;
			for (k=-1;k<2;k++)
			{
				for (l=-1;l<2;l++)
				{
					for (m=-1;m<2;m++)
					{
						pdx=dx+(k*ax)+(l*bx)+(m*cx);
						pdy=dy+(k*ay)+(l*by)+(m*cy);
						pdz=dz+(k*az)+(l*bz)+(m*cz);
						distance=sqrt((pdx*pdx)+(pdy*pdy)+(pdz*pdz));
						if (distance < *(double *)(dmatrix->data + i*dmatrix->strides[0]+j*dmatrix->strides[1]))
						{
							*(double *)(dmatrix->data + i*dmatrix->strides[0]+j*dmatrix->strides[1])=distance;
							*(double *)(dmatrix->data + j*dmatrix->strides[0]+i*dmatrix->strides[1])=distance;
						}
					}
				}
			}
		}
	}
	return PyArray_Return(dmatrix);
}


static PyObject *
geoext_blmatrix(PyObject *self, PyObject *args)
{
	PyArrayObject *blm;
	PyArrayObject *types;
	PyArrayObject *corads;
	PyObject *typelist, *coradlist;
	
	int natoms;
	int i,j;
	int dimensions[2];
	float rada,radb,bondlen;
	
	if (!PyArg_ParseTuple(args,"O!O!",&PyList_Type,&typelist,&PyList_Type,&coradlist))
		return NULL;
	
	natoms=PyList_Size(typelist);
	dimensions[0]=natoms;
	dimensions[1]=natoms;
	
	types =  (PyArrayObject *)PyArray_ContiguousFromObject(typelist,  PyArray_INT, 1,1);
	corads = (PyArrayObject *)PyArray_ContiguousFromObject(coradlist, PyArray_DOUBLE, 1,1);
	
	if (types==NULL || corads==NULL)
		return NULL;
	
	
	blm=(PyArrayObject *)PyArray_FromDims(2,dimensions,PyArray_DOUBLE);
	
	for (i=0;i<natoms;i++)
	{
		*(double *)(blm->data + i*blm->strides[0]+i*blm->strides[1])=-1000000.0;
		for (j=i+1;j<natoms;j++)
		{
			rada=*(double *)(corads->data + corads->strides[0] * 
				*(int *)(types->data + i * types->strides[0]));
			radb=*(double *)(corads->data + corads->strides[0] * 
				*(int *)(types->data + j * types->strides[0]));
			bondlen=rada+radb;
			*(double *)(blm->data + i*blm->strides[0]+j*blm->strides[1])=bondlen;
			*(double *)(blm->data + j*blm->strides[0]+i*blm->strides[1])=bondlen;
		}
	}
	return PyArray_Return(blm);
}



static PyObject *
geoext_blist(PyObject *self, PyObject *args)
{
	PyArrayObject *pos;
	PyListObject  *CORAD;
	PyListObject  *types;
	PyListObject  *bondlist;
	long natoms,typea,typeb;
	double rada,blen;
	long i,j;
	double dx,dy,dz,distance;
	double tolerance=1.1;
	
	if(!PyArg_ParseTuple(args, "O!O!O!|d", &PyArray_Type, &pos, &PyList_Type, &types, &PyList_Type, &CORAD, &tolerance))
		return NULL;
	if (pos->nd !=2 || pos->descr->type_num != PyArray_DOUBLE
			|| pos->dimensions[1]!=3){
		PyErr_SetString(PyExc_ValueError, "positions array must be two-dimensional, of shape (n,3) and of type float");
		return NULL;
	}
	natoms=pos->dimensions[0];
	/*initialize bond list as list of empty lists*/
	bondlist=PyList_New(natoms);
	for (i=0;i<natoms;i++)
	{
		PyList_SET_ITEM(	bondlist, i, PyList_New(0));
	}
	for (i=0;i<natoms;i++)
	{
		typea=PyInt_AsLong(PyList_GetItem(	types, i));
		rada=PyFloat_AsDouble(PyList_GetItem(	CORAD, typea));
		for (j=i+1;j<natoms;j++)
		{
			dx=*(double *)(pos->data + i * pos->strides[0])
				- *(double *)(pos->data + j * pos->strides[0]);
			dy=*(double *)(pos->data + i * pos->strides[0] + pos->strides[1])
				- *(double *)(pos->data + j * pos->strides[0] + pos->strides[1]);
			dz=*(double *)(pos->data + i * pos->strides[0] + 2 * pos->strides[1])
				- *(double *)(pos->data + j * pos->strides[0] + 2 * pos->strides[1]);
			distance=sqrt((dx*dx)+(dy*dy)+(dz*dz));
			typeb=PyInt_AsLong(PyList_GetItem(	types, j));
			blen=(tolerance* (double)PyFloat_AsDouble( PyList_GetItem(	CORAD, typeb))+rada)/5.291772e-01;
			if (distance < blen)
			{
				PyList_Append(PyList_GetItem(bondlist,i),PyInt_FromLong(j));
				PyList_Append(PyList_GetItem(bondlist,j),PyInt_FromLong(i));
			}
		}
	}
	return bondlist;
}


static PyObject *
geoext_sblist(PyObject *self, PyObject *args)
{
	PyArrayObject *pos;
	PyArrayObject *lattice;
	PyListObject  *CORAD;
	PyListObject  *types;
	PyListObject  *bondlist;
	long natoms,typea,typeb;
	double rada,blen;
	long i,j,k,l,m;
	double dx,dy,dz,distance;
	double pdx,pdy,pdz;
	double ax,ay,az,bx,by,bz,cx,cy,cz;
	double tolerance=1.1;
	
	
	if(!PyArg_ParseTuple(args, "O!O!O!O!|d", &PyArray_Type, &pos, &PyArray_Type,&lattice, &PyList_Type, &types, &PyList_Type, &CORAD, &tolerance))
		return NULL;
	if (pos->nd !=2 || pos->descr->type_num != PyArray_DOUBLE
			|| pos->dimensions[1]!=3){
		PyErr_SetString(PyExc_ValueError, "positions array must be two-dimensional, of shape (n,3) and of type float");
		return NULL;
	}
	if (lattice->nd !=2 || lattice->descr->type_num != PyArray_DOUBLE
			|| lattice->dimensions[0]!=3 || lattice->dimensions[1] !=3){
		PyErr_SetString(PyExc_ValueError, "lattice array must be two-dimensional, of shape (3,3) and of type float");
		return NULL;
	}
	natoms=pos->dimensions[0];
	/*initialize bond list as list of empty lists*/
	bondlist=PyList_New(natoms);
	for (i=0;i<natoms;i++)
	{
		PyList_SET_ITEM(	bondlist, i, PyList_New(0));
	}
	ax= *(double *)(lattice->data);
	ay= *(double *)(lattice->data+lattice->strides[1]);
	az= *(double *)(lattice->data+2*lattice->strides[1]);
	bx= *(double *)(lattice->data+lattice->strides[0]);
	by= *(double *)(lattice->data+lattice->strides[1]+lattice->strides[0]);
	bz= *(double *)(lattice->data+2*lattice->strides[1]+lattice->strides[0]);
	cx= *(double *)(lattice->data+2*lattice->strides[0]);
	cy= *(double *)(lattice->data+lattice->strides[1]+2*lattice->strides[0]);
	cz= *(double *)(lattice->data+2*lattice->strides[1]+2*lattice->strides[0]);
	for (i=0;i<natoms;i++)
	{
		typea=PyInt_AsLong(PyList_GetItem(	types, i));
		rada=PyFloat_AsDouble(PyList_GetItem(	CORAD, typea));
		for (j=i+1;j<natoms;j++)
		{
			dx=*(double *)(pos->data + i * pos->strides[0])
				- *(double *)(pos->data + j * pos->strides[0]);
			dy=*(double *)(pos->data + i * pos->strides[0] + pos->strides[1])
				- *(double *)(pos->data + j * pos->strides[0] + pos->strides[1]);
			dz=*(double *)(pos->data + i * pos->strides[0] + 2 * pos->strides[1])
				- *(double *)(pos->data + j * pos->strides[0] + 2 * pos->strides[1]);
			for (k=-1;k<2;k++)
			{
				for (l=-1;l<2;l++)
				{
					for (m=-1;m<2;m++)
					{
						pdx=dx+(k*ax)+(l*bx)+(m*cx);
						pdy=dy+(k*ay)+(l*by)+(m*cy);
						pdz=dz+(k*az)+(l*bz)+(m*cz);
						distance=sqrt((pdx*pdx)+(pdy*pdy)+(pdz*pdz));
						typeb=PyInt_AsLong(PyList_GetItem(	types, j));
						blen=(tolerance* (double)PyFloat_AsDouble( PyList_GetItem(	CORAD, typeb))+rada)/5.291772e-01;
						if (distance < blen)
						{
							PyList_Append(PyList_GetItem(bondlist,i),PyInt_FromLong(j));
							PyList_Append(PyList_GetItem(bondlist,j),PyInt_FromLong(i));
						}
					}
				}
			}
		}
	}
	return bondlist;
}


static PyMethodDef GeoExtMethods[] = {
	{"dmatrix", geoext_dmatrix, METH_VARARGS, 
		"return the (n,n) distance matrix for a numarray of (n,3) float coordinates"},
	{"sdmatrix", geoext_sdmatrix, METH_VARARGS, 
		"return the (n,n) minimum distance matrix for a numarray of (n,3) float coordinates with (3,3) array of lattice vectors"},
	{"blmatrix", geoext_blmatrix, METH_VARARGS,
		"return the (n,n) matrix of covalent radius sums for the types and covalent radius lists provided"},
	{"blist", geoext_blist, METH_VARARGS,
		"return the per-atom list of bond parter lists, molecule version"},
	{"sblist", geoext_sblist, METH_VARARGS,
		"return the per-atom list of bond parter lists, supercell version"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initgeoext(void)
{
    (void) Py_InitModule("geoext", GeoExtMethods);
	import_array();
}

