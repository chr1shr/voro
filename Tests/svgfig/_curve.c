#include <Python.h>

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

struct sample {
  double t, x, y;
  int discontinuity;
  struct sample *left;
  struct sample *right;
};

struct common_block {
  PyObject *parametric;
  PyObject *listoftrans;
  int counter;
  PyObject *random_sampling;
  int recursion_limit;
  double linearity_limit;
  double discontinuity_limit;
  unsigned int MT[624];
  int count624;
};

/* An implementation of the Mersenne Twistor random algorithm */
/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura. */
/* Copied from ROOT's TRandom3 */
/* (wanted to avoid compile-time dependencies for such a small part of the program) */
static void _curve_setseed(int seed, struct common_block *block) {
  block->count624 = 624;
  if (seed > 0) {
    block->MT[0] = seed;
  }
  int i;
  for (i = 1;  i < 624;  i++) {
    block->MT[i] = (1812433253 * (block->MT[i-1] ^ (block->MT[i-1] >> 30)) + i);
  }
}

static double _curve_random(struct common_block *block) {
   unsigned int y;
   int kM = 397;
   int kN = 624;
   unsigned int kTemperingMaskB = 0x9d2c5680;
   unsigned int kTemperingMaskC = 0xefc60000;
   unsigned int kUpperMask = 0x80000000;
   unsigned int kLowerMask = 0x7fffffff;
   unsigned int kMatrixA = 0x9908b0df;

   if (block->count624 >= kN) {
      int i;
      for (i = 0;  i < kN-kM;  i++) {
         y = (block->MT[i] & kUpperMask) | (block->MT[i+1] & kLowerMask);
         block->MT[i] = block->MT[i+kM] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      }

      for (;  i < kN-1;  i++) {
         y = (block->MT[i] & kUpperMask) | (block->MT[i+1] & kLowerMask);
         block->MT[i] = block->MT[i+kM-kN] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      }

      y = (block->MT[kN-1] & kUpperMask) | (block->MT[0] & kLowerMask);
      block->MT[kN-1] = block->MT[kM-1] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      block->count624 = 0;
   }

   y = block->MT[block->count624++];
   y ^= (y >> 11);
   y ^= ((y << 7 ) & kTemperingMaskB);
   y ^= ((y << 15) & kTemperingMaskC);
   y ^= (y >> 18);

   if (y != 0) return ((double) y * 2.3283064365386963e-10); // * pow(2, -32)
   return _curve_random(block);
}

/* evaluate parametric function at a point, passing it through a list of coordinate transformations */
static int _curve_eval(PyObject *parametric, PyObject *listoftrans, double t, double *fx, double *fy) {
  PyObject *args = Py_BuildValue("(d)", t);
  PyObject *result = PyObject_CallObject(parametric, args);
  if (result == NULL) {
    Py_DECREF(args);
    return 0;
  }
  Py_DECREF(args);

  if (!PySequence_Check(result)  ||  PySequence_Size(result) != 2) {
    PyErr_SetString(PyExc_TypeError, "The parametric function must return two real values.");
    Py_DECREF(result);
    return 0;
  }
  
  PyObject *x = PySequence_GetItem(result, 0);
  PyObject *y = PySequence_GetItem(result, 1);
  Py_DECREF(result);

  if (!PyNumber_Check(x)  ||  !PyNumber_Check(y)) {
    PyErr_SetString(PyExc_TypeError, "The parametric function must return two real values.");
    Py_DECREF(x);
    Py_DECREF(y);
    return 0;
  }

  int i;
  int lenlistoftrans = PySequence_Size(listoftrans);
  for (i = 0;  i < lenlistoftrans;  i++) {
    PyObject *trans = PySequence_GetItem(listoftrans, i);
    
    args = Py_BuildValue("(OO)", x, y);
    Py_DECREF(x);
    Py_DECREF(y);
    result = PyObject_CallObject(trans, args);
    if (result == NULL) {
      Py_DECREF(trans);
      Py_DECREF(args);
      return 0;
    }
    Py_DECREF(trans);
    Py_DECREF(args);

    if (!PySequence_Check(result)  ||  PySequence_Size(result) != 2) {
      PyErr_SetString(PyExc_TypeError, "The transformation functions must return two real values.");
      Py_DECREF(result);
      return 0;
    }

    x = PySequence_GetItem(result, 0);  // only borrowed
    y = PySequence_GetItem(result, 1);  // only borrowed
    Py_DECREF(result);

    if (!PyNumber_Check(x)  ||  !PyNumber_Check(y)) {
      PyErr_SetString(PyExc_TypeError, "The transformation functions must return two real values.");
      Py_DECREF(x);
      Py_DECREF(y);
      return 0;
    }
  }

  *fx = PyFloat_AsDouble(x);
  *fy = PyFloat_AsDouble(y);
  Py_DECREF(x);
  Py_DECREF(y);

  return 1;
}

/* recursively called to fill a (doubly-linked) list of sample points where it needs it most */
/* with the default parameters, it computes a few more points than are typically needed */
/* after pruning the extras, that guarantees a nice smooth curve */
int _curve_subsample(struct sample *left, struct sample *right, int depth, struct common_block *block) {
  /* make new mid node and link it up */
  struct sample *mid = (struct sample*)malloc(sizeof(struct sample));
  if (mid == NULL) {
    PyErr_Format(PyExc_MemoryError, "Ran out of memory while sampling function (%d nodes created)", block->counter);
    return 0;
  }
  block->counter++;

  if (block->random_sampling == Py_True) {
    mid->t = left->t + (0.3 + 0.4*_curve_random(block))*(right->t - left->t);
  }
  else {
    mid->t = left->t + 0.5*(right->t - left->t);
  }
  
  if (!_curve_eval(block->parametric, block->listoftrans, mid->t, &(mid->x), &(mid->y))) {
    free(mid);
    return 0;
  }
  mid->discontinuity = 0;

  left->right = mid;
  mid->left = left;
  mid->right = right;
  right->left = mid;

  /* calculate the distance of closest approach of mid to the line between left and right */
  double numer = (left->x)*(right->y - mid->y) + (mid->x)*(left->y - right->y) + (right->x)*(mid->y - left->y);
  double denom = sqrt((left->x - right->x)*(left->x - right->x) + (left->y - right->y)*(left->y - right->y));

  if (depth < 3  ||
      (denom == 0.  &&  left->t != right->t)  ||
      (denom > block->discontinuity_limit)  ||
      (denom != 0.  &&  fabs(numer/denom) > block->linearity_limit)) {

    if (depth < block->recursion_limit) {
      if (!_curve_subsample(left, mid, depth+1, block)) return 0;
      if (!_curve_subsample(mid, right, depth+1, block)) return 0;
    }

    else {
      /* we've sampled many points and yet it's still not a small linear gap */
      /* break the line: we've found a discontinuity */
      mid->discontinuity = 1;
    }
  }

  return 1;
}

/* the only function which is called from the outside: the interface to Python */
static PyObject *_curve_curve(PyObject *self, PyObject *args, PyObject *kwds) {
  const char *errstring = "arguments are: parametric function to plot, list of transformations to apply to each point, low endpoint, high endpoint.  \nkeyword arguments are: random_sampling (True), random_seed (12345), recursion_limit (15), linearity_limit (0.05), discontinuity_limit (5.)";

  PyObject *parametric;
  PyObject *listoftrans;
  double low, high;
  PyObject *random_sampling = Py_True;
  int random_seed = 12345;
  int recursion_limit = 15;
  double linearity_limit = 0.05;
  double discontinuity_limit = 5.;

  static char *kwlist[] = {"parametric", "listoftrans", "low", "high", "random_sampling", "random_seed", "recursion_limit", "linearity_limit", "discontinuity_limit", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOdd|Oiidd", kwlist, &parametric, &listoftrans, &low, &high, &random_sampling, &random_seed, &recursion_limit, &linearity_limit, &discontinuity_limit)) {
    PyErr_SetString(PyExc_TypeError, errstring);
    return NULL;
  }

  if (!PyCallable_Check(parametric)) {
    PyErr_SetString(PyExc_TypeError, errstring);
    return NULL;
  }

  if (!PySequence_Check(listoftrans)) {
    PyErr_SetString(PyExc_TypeError, errstring);
    return NULL;
  }

  int i;
  int lenlistoftrans = PySequence_Size(listoftrans);
  for (i = 0;  i < lenlistoftrans;  i++) {
    PyObject *trans = PySequence_GetItem(listoftrans, i);
    if (!PyCallable_Check(trans)) {
      PyErr_SetString(PyExc_TypeError, errstring);
      Py_DECREF(trans);
      return NULL;
    }
    Py_DECREF(trans);
  }

  if (random_sampling != Py_True  &&  random_sampling != Py_False) {
    PyErr_SetString(PyExc_TypeError, errstring);
    return NULL;
  }

  /* build doubly-linked list of samples */
  struct sample *samplelow = (struct sample*)malloc(sizeof(struct sample));
  if (samplelow == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Ran out of memory while sampling function (1 node created)");
    return NULL;
  }
  samplelow->t = low;
  if (!_curve_eval(parametric, listoftrans, samplelow->t, &(samplelow->x), &(samplelow->y))) {
    free(samplelow);
    return NULL;
  }
  samplelow->discontinuity = 0;

  struct sample *samplehigh = (struct sample*)malloc(sizeof(struct sample));
  if (samplelow == NULL) {
    PyErr_SetString(PyExc_MemoryError, "Ran out of memory while sampling function (2 nodes created)");
    free(samplelow);
    return NULL;
  }
  samplehigh->t = high;
  if (!_curve_eval(parametric, listoftrans, samplehigh->t, &(samplehigh->x), &(samplehigh->y))) {
    free(samplelow);
    free(samplehigh);
    return NULL;
  }
  samplehigh->discontinuity = 0;

  samplelow->left = NULL;
  samplelow->right = samplehigh;
  samplehigh->left = samplelow;
  samplehigh->right = NULL;

  struct common_block block;
  block.parametric = parametric;
  block.listoftrans = listoftrans;
  block.counter = 2;
  block.random_sampling = random_sampling;
  block.recursion_limit = recursion_limit;
  block.linearity_limit = linearity_limit;
  block.discontinuity_limit = discontinuity_limit;
  _curve_setseed(random_seed, &block);

  /* recursively find most of the points */
  if (!_curve_subsample(samplelow, samplehigh, 0, &block)) {
    /* free nodes in case of error */
    struct sample *p = samplelow;
    while (p != NULL) {
      struct sample *next = p->right;
      free(p);
      p = next;
    }
    return NULL;
  }
  
  /* prune excess points that are within the linearity bounds */
  struct sample *left = samplelow;
  int length = 1; /* get the post-pruning length */
  while (left->right != NULL) {
    struct sample *mid = left->right;
    struct sample *right = mid->right;

    if (right != NULL  &&  !left->discontinuity  &&  !mid->discontinuity  &&  !right->discontinuity) {
      double numer = (left->x)*(right->y - mid->y) + (mid->x)*(left->y - right->y) + (right->x)*(mid->y - left->y);
      double denom = sqrt((left->x - right->x)*(left->x - right->x) + (left->y - right->y)*(left->y - right->y));

      if (denom != 0.  &&  fabs(numer/denom) < linearity_limit) {
	free(mid); /* drop this point; it doesn't contribute to the smoothness of the curve */
	left->right = right;
	right->left = left;
      }
      else {
	length++;
	left = left->right; /* increment left */
      }
    }
    else {
      length++;
      left = left->right; /* increment left */
    }
  }

  /* return a Python tuple of numbers and free all nodes */
  PyObject *output = PyTuple_New(length);
  int failure = 0;

  struct sample *p = samplelow;
  for (i = 0;  i < length;  i++) {
    if (p->discontinuity) {
      Py_INCREF(Py_None);
      if (PyTuple_SetItem(output, i, Py_None) != 0) failure = 1;
    }
    else {
      if (PyTuple_SetItem(output, i, Py_BuildValue("dd", p->x, p->y)) != 0) failure = 1;
    }

    struct sample *next = p->right;
    free(p);
    p = next;
  }

  /* it's important to finish freeing nodes first */
  if (failure) {
    Py_DECREF(output);
    return NULL;
  }

  return output;
}

static PyMethodDef _curve_methods[] = {
  {"curve", ((PyCFunction)(_curve_curve)), METH_VARARGS | METH_KEYWORDS, ""},
  {NULL}
};

PyMODINIT_FUNC init_curve() {
  Py_InitModule3("_curve", _curve_methods, "");
}
