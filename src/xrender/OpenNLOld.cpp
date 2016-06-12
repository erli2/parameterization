/*
 *  $Id: opennlOLD.c,v 1.3 2005/10/30 18:38:35 blendix Exp $
 *
 *  OpenNLOLD: Numerical Library
 *  Copyright (C) 2004 Bruno Levy
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ISA Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#include "stdafx.h"

#include "OpenNLOld.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef NLOLD_PARANOID
#ifndef NLOLD_DEBUG
#define NLOLD_DEBUG
#endif
#endif

/* SuperLU includes */
#include "superlu/ssp_defs.h"
#include "superlu/util.h"

/************************************************************************************/
/* Assertions */


static void __nlOLD_assertion_failed(char* cond, char* file, int line) {
    fprintf(
        stderr, 
        "OpenNLOLD assertion failed: %s, file:%s, line:%d\n",
        cond,file,line
    );
    abort();
}

static void __nlOLD_range_assertion_failed(
    float x, float min_val, float max_val, char* file, int line
) {
    fprintf(
        stderr, 
        "OpenNLOLD range assertion failed: %f in [ %f ... %f ], file:%s, line:%d\n",
        x, min_val, max_val, file,line
    );
    abort();
}

static void __nlOLD_should_not_have_reached(char* file, int line) {
    fprintf(
        stderr, 
        "OpenNLOLD should not have reached this point: file:%s, line:%d\n",
        file,line
    );
    abort();
}


#define __nlOLD_assert(x) {                                        \
    if(!(x)) {                                                  \
        __nlOLD_assertion_failed(#x,__FILE__, __LINE__);          \
    }                                                           \
} 

#define __nlOLD_range_assert(x,min_val,max_val) {                  \
    if(((x) < (min_val)) || ((x) > (max_val))) {                \
        __nlOLD_range_assertion_failed(x, min_val, max_val,        \
            __FILE__, __LINE__                                  \
        );                                                     \
    }                                                           \
}

#define __nlOLD_assert_not_reached {                               \
    __nlOLD_should_not_have_reached(__FILE__, __LINE__);          \
}

//#ifdef NLOLD_DEBUG
//#define __nlOLD_debug_assert(x) __nlOLD_assert(x)
//#define __nlOLD_debug_range_assert(x,min_val,max_val) __nlOLD_range_assert(x,min_val,max_val)
//#else
#define __nlOLD_debug_assert(x) 
#define __nlOLD_debug_range_assert(x,min_val,max_val) 
//#endif

//#ifdef NLOLD_PARANOID
//#define __nlOLD_parano_assert(x) __nlOLD_assert(x)
//#define __nlOLD_parano_range_assert(x,min_val,max_val) __nlOLD_range_assert(x,min_val,max_val)
//#else
#define __nlOLD_parano_assert(x) 
#define __nlOLD_parano_range_assert(x,min_val,max_val) 
//#endif

/************************************************************************************/
/* classic macros */

#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y)) 
#endif

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y)) 
#endif

/************************************************************************************/
/* memory management */

#define __NLOLD_NEW(T)                (T*)(calloc(1, sizeof(T))) 
#define __NLOLD_NEW_ARRAY(T,NB)       (T*)(calloc((NB),sizeof(T))) 
#define __NLOLD_RENEW_ARRAY(T,x,NB)   (T*)(realloc(x,(NB)*sizeof(T))) 
#define __NLOLD_DELETE(x)             free(x); x = NULL 
#define __NLOLD_DELETE_ARRAY(x)       free(x); x = NULL

#define __NLOLD_CLEAR(T, x)           memset(x, 0, sizeof(T)) 
#define __NLOLD_CLEAR_ARRAY(T,x,NB)   memset(x, 0, (NB)*sizeof(T)) 

/************************************************************************************/
/* Dynamic arrays for sparse row/columns */

typedef struct {
    NLOLDuint   index;
    NLOLDfloat value;
} __NLOLDCoeff;

typedef struct {
    NLOLDuint size;
    NLOLDuint capacity;
    __NLOLDCoeff* coeff;
} __NLOLDRowColumn;

static void __nlOLDRowColumnConstruct(__NLOLDRowColumn* c) {
    c->size     = 0;
    c->capacity = 0;
    c->coeff    = NULL;
}

static void __nlOLDRowColumnDestroy(__NLOLDRowColumn* c) {
    __NLOLD_DELETE_ARRAY(c->coeff);
#ifdef NLOLD_PARANOID
    __NLOLD_CLEAR(__NLOLDRowColumn, c); 
#endif
}

static void __nlOLDRowColumnGrow(__NLOLDRowColumn* c) {
    if(c->capacity != 0) {
        c->capacity = 2 * c->capacity;
        c->coeff = __NLOLD_RENEW_ARRAY(__NLOLDCoeff, c->coeff, c->capacity);
    } else {
        c->capacity = 4;
        c->coeff = __NLOLD_NEW_ARRAY(__NLOLDCoeff, c->capacity);
    }
}

static void __nlOLDRowColumnAdd(__NLOLDRowColumn* c, NLOLDint index, NLOLDfloat value) {
    NLOLDuint i;
    for(i=0; i<c->size; i++) {
        if(c->coeff[i].index == (NLOLDuint)index) {
            c->coeff[i].value += value;
            return;
        }
    }
    if(c->size == c->capacity) {
        __nlOLDRowColumnGrow(c);
    }
    c->coeff[c->size].index = index;
    c->coeff[c->size].value = value;
    c->size++;
}

/* Does not check whether the index already exists */
static void __nlOLDRowColumnAppend(__NLOLDRowColumn* c, NLOLDint index, NLOLDfloat value) {
    if(c->size == c->capacity) {
        __nlOLDRowColumnGrow(c);
    }
    c->coeff[c->size].index = index;
    c->coeff[c->size].value = value;
    c->size++;
}

static void __nlOLDRowColumnZero(__NLOLDRowColumn* c) {
    c->size = 0;
}

static void __nlOLDRowColumnClear(__NLOLDRowColumn* c) {
    c->size     = 0;
    c->capacity = 0;
    __NLOLD_DELETE_ARRAY(c->coeff);
}

/************************************************************************************/
/* SparseMatrix data structure */

#define __NLOLD_ROWS      1
#define __NLOLD_COLUMNS   2
#define __NLOLD_SYMMETRIC 4

typedef struct {
    NLOLDuint m;
    NLOLDuint n;
    NLOLDuint diag_size;
    NLOLDenum storage;
    __NLOLDRowColumn* row;
    __NLOLDRowColumn* column;
    NLOLDfloat*      diag;
} __NLOLDSparseMatrix;


static void __nlOLDSparseMatrixConstruct(
    __NLOLDSparseMatrix* M, NLOLDuint m, NLOLDuint n, NLOLDenum storage
) {
    NLOLDuint i;
    M->m = m;
    M->n = n;
    M->storage = storage;
    if(storage & __NLOLD_ROWS) {
        M->row = __NLOLD_NEW_ARRAY(__NLOLDRowColumn, m);
        for(i=0; i<n; i++) {
            __nlOLDRowColumnConstruct(&(M->row[i]));
        }
    } else {
        M->row = NULL;
    }

    if(storage & __NLOLD_COLUMNS) {
        M->column = __NLOLD_NEW_ARRAY(__NLOLDRowColumn, n);
        for(i=0; i<n; i++) {
            __nlOLDRowColumnConstruct(&(M->column[i]));
        }
    } else {
        M->column = NULL;
    }

    M->diag_size = MIN(m,n);
    M->diag = __NLOLD_NEW_ARRAY(NLOLDfloat, M->diag_size);
}

static void __nlOLDSparseMatrixDestroy(__NLOLDSparseMatrix* M) {
    NLOLDuint i;
    __NLOLD_DELETE_ARRAY(M->diag);
    if(M->storage & __NLOLD_ROWS) {
        for(i=0; i<M->m; i++) {
            __nlOLDRowColumnDestroy(&(M->row[i]));
        }
        __NLOLD_DELETE_ARRAY(M->row);
    }
    if(M->storage & __NLOLD_COLUMNS) {
        for(i=0; i<M->n; i++) {
            __nlOLDRowColumnDestroy(&(M->column[i]));
        }
        __NLOLD_DELETE_ARRAY(M->column);
    }
#ifdef NLOLD_PARANOID
    __NLOLD_CLEAR(__NLOLDSparseMatrix,M);
#endif
}

static void __nlOLDSparseMatrixAdd(
    __NLOLDSparseMatrix* M, NLOLDuint i, NLOLDuint j, NLOLDfloat value
) {
    __nlOLD_parano_range_assert(i, 0, M->m - 1);
    __nlOLD_parano_range_assert(j, 0, M->n - 1);
    if((M->storage & __NLOLD_SYMMETRIC) && (j > i)) {
        return;
    }
    if(i == j) {
        M->diag[i] += value;
    }
    if(M->storage & __NLOLD_ROWS) {
        __nlOLDRowColumnAdd(&(M->row[i]), j, value);
    }
    if(M->storage & __NLOLD_COLUMNS) {
        __nlOLDRowColumnAdd(&(M->column[j]), i, value);
    }
}

static void __nlOLDSparseMatrixClear( __NLOLDSparseMatrix* M) {
    NLOLDuint i;
    if(M->storage & __NLOLD_ROWS) {
        for(i=0; i<M->m; i++) {
            __nlOLDRowColumnClear(&(M->row[i]));
        }
    }
    if(M->storage & __NLOLD_COLUMNS) {
        for(i=0; i<M->n; i++) {
            __nlOLDRowColumnClear(&(M->column[i]));
        }
    }
    __NLOLD_CLEAR_ARRAY(NLOLDfloat, M->diag, M->diag_size);    
}

/* Returns the number of non-zero coefficients */
static NLOLDuint __nlOLDSparseMatrixNNZ( __NLOLDSparseMatrix* M) {
    NLOLDuint nnz = 0;
    NLOLDuint i;
    if(M->storage & __NLOLD_ROWS) {
        for(i = 0; i<M->m; i++) {
            nnz += M->row[i].size;
        }
    } else if (M->storage & __NLOLD_COLUMNS) {
        for(i = 0; i<M->n; i++) {
            nnz += M->column[i].size;
        }
    } else {
        __nlOLD_assert_not_reached;
    }
    return nnz;
}

/************************************************************************************/
/* SparseMatrix x Vector routines, internal helper routines */

static void __nlOLDSparseMatrix_mult_rows_symmetric(
    __NLOLDSparseMatrix* A, NLOLDfloat* x, NLOLDfloat* y
) {
    NLOLDuint m = A->m;
    NLOLDuint i,ij;
    __NLOLDRowColumn* Ri = NULL;
    __NLOLDCoeff* c = NULL;
    for(i=0; i<m; i++) {
        y[i] = 0;
        Ri = &(A->row[i]);
        for(ij=0; ij<Ri->size; ij++) {
            c = &(Ri->coeff[ij]);
            y[i] += c->value * x[c->index];
            if(i != c->index) {
                y[c->index] += c->value * x[i];
            }
        }
    }
}

static void __nlOLDSparseMatrix_mult_rows(
    __NLOLDSparseMatrix* A, NLOLDfloat* x, NLOLDfloat* y
) {
    NLOLDuint m = A->m;
    NLOLDuint i,ij;
    __NLOLDRowColumn* Ri = NULL;
    __NLOLDCoeff* c = NULL;
    for(i=0; i<m; i++) {
        y[i] = 0;
        Ri = &(A->row[i]);
        for(ij=0; ij<Ri->size; ij++) {
            c = &(Ri->coeff[ij]);
            y[i] += c->value * x[c->index];
        }
    }
}

static void __nlOLDSparseMatrix_mult_cols_symmetric(
    __NLOLDSparseMatrix* A, NLOLDfloat* x, NLOLDfloat* y
) {
    NLOLDuint n = A->n;
    NLOLDuint j,ii;
    __NLOLDRowColumn* Cj = NULL;
    __NLOLDCoeff* c = NULL;
    for(j=0; j<n; j++) {
        y[j] = 0;
        Cj = &(A->column[j]);
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]);
            y[c->index] += c->value * x[j];
            if(j != c->index) {
                y[j] += c->value * x[c->index];
            }
        }
    }
}

static void __nlOLDSparseMatrix_mult_cols(
    __NLOLDSparseMatrix* A, NLOLDfloat* x, NLOLDfloat* y
) {
    NLOLDuint n = A->n;
    NLOLDuint j,ii; 
    __NLOLDRowColumn* Cj = NULL;
    __NLOLDCoeff* c = NULL;
    __NLOLD_CLEAR_ARRAY(NLOLDfloat, y, A->m);
    for(j=0; j<n; j++) {
        Cj = &(A->column[j]);
        for(ii=0; ii<Cj->size; ii++) {
            c = &(Cj->coeff[ii]);
            y[c->index] += c->value * x[j];
        }
    }
}

/************************************************************************************/
/* SparseMatrix x Vector routines, main driver routine */

static void __nlOLDSparseMatrixMult(__NLOLDSparseMatrix* A, NLOLDfloat* x, NLOLDfloat* y) {
    if(A->storage & __NLOLD_ROWS) {
        if(A->storage & __NLOLD_SYMMETRIC) {
            __nlOLDSparseMatrix_mult_rows_symmetric(A, x, y);
        } else {
            __nlOLDSparseMatrix_mult_rows(A, x, y);
        }
    } else {
        if(A->storage & __NLOLD_SYMMETRIC) {
            __nlOLDSparseMatrix_mult_cols_symmetric(A, x, y);
        } else {
            __nlOLDSparseMatrix_mult_cols(A, x, y);
        }
    }
}

/************************************************************************************/
/* NLOLDContext data structure */

typedef void(*__NLOLDMatrixFunc)(float* x, float* y);

typedef struct {
    NLOLDfloat  value;
    NLOLDboolean locked;
    NLOLDuint    index;
} __NLOLDVariable;

#define __NLOLD_STATE_INITIAL            0
#define __NLOLD_STATE_SYSTEM             1
#define __NLOLD_STATE_MATRIX             2
#define __NLOLD_STATE_ROW                3
#define __NLOLD_STATE_MATRIX_CONSTRUCTED 4
#define __NLOLD_STATE_SYSTEM_CONSTRUCTED 5
#define __NLOLD_STATE_SOLVED             6

typedef struct {
    NLOLDenum           state;
    __NLOLDVariable*    variable;
    NLOLDuint           n;
    __NLOLDSparseMatrix M;
    __NLOLDRowColumn    af;
    __NLOLDRowColumn    al;
    __NLOLDRowColumn    xl;
    NLOLDfloat*        x;
    NLOLDfloat*        b;
    NLOLDfloat         right_hand_side;
    NLOLDfloat         row_scaling;
    NLOLDuint           nb_variables;
    NLOLDuint           current_row;
    NLOLDboolean        least_squares;
    NLOLDboolean        symmetric;
    NLOLDboolean        normalize_rows;
    NLOLDboolean        alloc_M;
    NLOLDboolean        alloc_af;
    NLOLDboolean        alloc_al;
    NLOLDboolean        alloc_xl;
    NLOLDboolean        alloc_variable;
    NLOLDboolean        alloc_x;
    NLOLDboolean        alloc_b;
    NLOLDfloat         error;
    __NLOLDMatrixFunc   matrix_vector_prod;
} __NLOLDContext;

static __NLOLDContext* __nlOLDCurrentContext = NULL;

static void __nlOLDMatrixVectorProd_default(NLOLDfloat* x, NLOLDfloat* y) {
    __nlOLDSparseMatrixMult(&(__nlOLDCurrentContext->M), x, y);
}


NLOLDContext nlOLDNewContext(void) {
    __NLOLDContext* result      = __NLOLD_NEW(__NLOLDContext);
    result->state            = __NLOLD_STATE_INITIAL;
    result->row_scaling      = 1.0;
    result->right_hand_side  = 0.0;
    result->matrix_vector_prod = __nlOLDMatrixVectorProd_default;
    nlOLDMakeCurrent(result);
    return result;
}

void nlOLDDeleteContext(NLOLDContext context_in) {
    __NLOLDContext* context = (__NLOLDContext*)(context_in);
    if(__nlOLDCurrentContext == context) {
        __nlOLDCurrentContext = NULL;
    }
    if(context->alloc_M) {
        __nlOLDSparseMatrixDestroy(&context->M);
    }
    if(context->alloc_af) {
        __nlOLDRowColumnDestroy(&context->af);
    }
    if(context->alloc_al) {
        __nlOLDRowColumnDestroy(&context->al);
    }
    if(context->alloc_xl) {
        __nlOLDRowColumnDestroy(&context->xl);
    }
    if(context->alloc_variable) {
        __NLOLD_DELETE_ARRAY(context->variable);
    }
    if(context->alloc_x) {
        __NLOLD_DELETE_ARRAY(context->x);
    }
    if(context->alloc_b) {
        __NLOLD_DELETE_ARRAY(context->b);
    }

#ifdef NLOLD_PARANOID
    __NLOLD_CLEAR(__NLOLDContext, context);
#endif
    __NLOLD_DELETE(context);
}

void nlOLDMakeCurrent(NLOLDContext context) {
    __nlOLDCurrentContext = (__NLOLDContext*)(context);
}

NLOLDContext nlOLDGetCurrent(void) {
    return __nlOLDCurrentContext;
}

static void __nlOLDCheckState(NLOLDenum state) {
    __nlOLD_assert(__nlOLDCurrentContext->state == state);
}

static void __nlOLDTransition(NLOLDenum from_state, NLOLDenum to_state) {
    __nlOLDCheckState(from_state);
    __nlOLDCurrentContext->state = to_state;
}

/************************************************************************************/
/* Get/Set parameters */

void nlOLDSolverParameterf(NLOLDenum pname, NLOLDfloat param) {
    __nlOLDCheckState(__NLOLD_STATE_INITIAL);
    switch(pname) {
    case NLOLD_NB_VARIABLES: {
        __nlOLD_assert(param > 0);
        __nlOLDCurrentContext->nb_variables = (NLOLDuint)param;
    } break;
    case NLOLD_LEAST_SQUARES: {
        __nlOLDCurrentContext->least_squares = (NLOLDboolean)param;
    } break;
    case NLOLD_SYMMETRIC: {
        __nlOLDCurrentContext->symmetric = (NLOLDboolean)param;        
    }
    default: {
        __nlOLD_assert_not_reached;
    } break;
    }
}

void nlOLDSolverParameteri(NLOLDenum pname, NLOLDint param) {
    __nlOLDCheckState(__NLOLD_STATE_INITIAL);
    switch(pname) {
    case NLOLD_NB_VARIABLES: {
        __nlOLD_assert(param > 0);
        __nlOLDCurrentContext->nb_variables = (NLOLDuint)param;
    } break;
    case NLOLD_LEAST_SQUARES: {
        __nlOLDCurrentContext->least_squares = (NLOLDboolean)param;
    } break;
    case NLOLD_SYMMETRIC: {
        __nlOLDCurrentContext->symmetric = (NLOLDboolean)param;        
    }
    default: {
        __nlOLD_assert_not_reached;
    } break;
    }
}

void nlOLDRowParameterf(NLOLDenum pname, NLOLDfloat param) {
    __nlOLDCheckState(__NLOLD_STATE_MATRIX);
    switch(pname) {
    case NLOLD_RIGHT_HAND_SIDE: {
        __nlOLDCurrentContext->right_hand_side = param;
    } break;
    case NLOLD_ROW_SCALING: {
        __nlOLDCurrentContext->row_scaling = param;
    } break;
    }
}

void nlOLDRowParameteri(NLOLDenum pname, NLOLDint param) {
    __nlOLDCheckState(__NLOLD_STATE_MATRIX);
    switch(pname) {
    case NLOLD_RIGHT_HAND_SIDE: {
        __nlOLDCurrentContext->right_hand_side = (NLOLDfloat)param;
    } break;
    case NLOLD_ROW_SCALING: {
        __nlOLDCurrentContext->row_scaling = (NLOLDfloat)param;
    } break;
    }
}

void nlOLDGetBooleanv(NLOLDenum pname, NLOLDboolean* params) {
    switch(pname) {
    case NLOLD_LEAST_SQUARES: {
        *params = __nlOLDCurrentContext->least_squares;
    } break;
    case NLOLD_SYMMETRIC: {
        *params = __nlOLDCurrentContext->symmetric;
    } break;
    default: {
        __nlOLD_assert_not_reached;
    } break;
    }
}

void nlOLDGetFloatv(NLOLDenum pname, NLOLDfloat* params) {
    switch(pname) {
    case NLOLD_NB_VARIABLES: {
        *params = (NLOLDfloat)(__nlOLDCurrentContext->nb_variables);
    } break;
    case NLOLD_LEAST_SQUARES: {
        *params = (NLOLDfloat)(__nlOLDCurrentContext->least_squares);
    } break;
    case NLOLD_SYMMETRIC: {
        *params = (NLOLDfloat)(__nlOLDCurrentContext->symmetric);
    } break;
    case NLOLD_ERROR: {
        *params = (NLOLDfloat)(__nlOLDCurrentContext->error);
    } break;
    default: {
        __nlOLD_assert_not_reached;
    } break;
    }
}

void nlOLDGetIntergerv(NLOLDenum pname, NLOLDint* params) {
    switch(pname) {
    case NLOLD_NB_VARIABLES: {
        *params = (NLOLDint)(__nlOLDCurrentContext->nb_variables);
    } break;
    case NLOLD_LEAST_SQUARES: {
        *params = (NLOLDint)(__nlOLDCurrentContext->least_squares);
    } break;
    case NLOLD_SYMMETRIC: {
        *params = (NLOLDint)(__nlOLDCurrentContext->symmetric);
    } break;
    default: {
        __nlOLD_assert_not_reached;
    } break;
    }
}

/************************************************************************************/
/* Enable / Disable */

void nlOLDEnable(NLOLDenum pname) {
    switch(pname) {
    case NLOLD_NORMALIZE_ROWS: {
        __nlOLD_assert(__nlOLDCurrentContext->state != __NLOLD_STATE_ROW);
        __nlOLDCurrentContext->normalize_rows = NLOLD_TRUE;
    } break;
    default: {
        __nlOLD_assert_not_reached;
    }
    }
}

void nlOLDDisable(NLOLDenum pname) {
    switch(pname) {
    case NLOLD_NORMALIZE_ROWS: {
        __nlOLD_assert(__nlOLDCurrentContext->state != __NLOLD_STATE_ROW);
        __nlOLDCurrentContext->normalize_rows = NLOLD_FALSE;
    } break;
    default: {
        __nlOLD_assert_not_reached;
    }
    }
}

NLOLDboolean nlOLDIsEnabled(NLOLDenum pname) {
    switch(pname) {
    case NLOLD_NORMALIZE_ROWS: {
        return __nlOLDCurrentContext->normalize_rows;
    } break;
    default: {
        __nlOLD_assert_not_reached;
    }
    }
    return NLOLD_FALSE;
}

/************************************************************************************/
/* Get/Set Lock/UnlOLDock variables */

void nlOLDSetVariable(NLOLDuint index, NLOLDfloat value) {
    __nlOLDCheckState(__NLOLD_STATE_SYSTEM);
    __nlOLD_parano_range_assert(index, 0, __nlOLDCurrentContext->nb_variables - 1);
    __nlOLDCurrentContext->variable[index].value = value;    
}

NLOLDfloat nlOLDGetVariable(NLOLDuint index) {
    __nlOLD_assert(__nlOLDCurrentContext->state != __NLOLD_STATE_INITIAL);
    __nlOLD_parano_range_assert(index, 0, __nlOLDCurrentContext->nb_variables - 1);
    return __nlOLDCurrentContext->variable[index].value;
}

void nlOLDLockVariable(NLOLDuint index) {
    __nlOLDCheckState(__NLOLD_STATE_SYSTEM);
    __nlOLD_parano_range_assert(index, 0, __nlOLDCurrentContext->nb_variables - 1);
    __nlOLDCurrentContext->variable[index].locked = NLOLD_TRUE;
}

void nlOLDUnlOLDockVariable(NLOLDuint index) {
    __nlOLDCheckState(__NLOLD_STATE_SYSTEM);
    __nlOLD_parano_range_assert(index, 0, __nlOLDCurrentContext->nb_variables - 1);
    __nlOLDCurrentContext->variable[index].locked = NLOLD_FALSE;
}

NLOLDboolean nlOLDVariableIsLocked(NLOLDuint index) {
    __nlOLD_assert(__nlOLDCurrentContext->state != __NLOLD_STATE_INITIAL);
    __nlOLD_parano_range_assert(index, 0, __nlOLDCurrentContext->nb_variables - 1);
    return __nlOLDCurrentContext->variable[index].locked;
}

/************************************************************************************/
/* System construction */

static void __nlOLDVariablesToVector() {
    NLOLDuint i;
    __nlOLD_assert(__nlOLDCurrentContext->alloc_x);
    __nlOLD_assert(__nlOLDCurrentContext->alloc_variable);
    for(i=0; i<__nlOLDCurrentContext->nb_variables; i++) {
        __NLOLDVariable* v = &(__nlOLDCurrentContext->variable[i]);
        if(!v->locked) {
            __nlOLD_assert(v->index < __nlOLDCurrentContext->n);
            __nlOLDCurrentContext->x[v->index] = v->value;
        }
    }
}

static void __nlOLDVectorToVariables() {
    NLOLDuint i;
    __nlOLD_assert(__nlOLDCurrentContext->alloc_x);
    __nlOLD_assert(__nlOLDCurrentContext->alloc_variable);
    for(i=0; i<__nlOLDCurrentContext->nb_variables; i++) {
        __NLOLDVariable* v = &(__nlOLDCurrentContext->variable[i]);
        if(!v->locked) {
            __nlOLD_assert(v->index < __nlOLDCurrentContext->n);
            v->value = __nlOLDCurrentContext->x[v->index];
        }
    }
}


static void __nlOLDBeginSystem() {
    __nlOLDTransition(__NLOLD_STATE_INITIAL, __NLOLD_STATE_SYSTEM);
    __nlOLD_assert(__nlOLDCurrentContext->nb_variables > 0);
    __nlOLDCurrentContext->variable = __NLOLD_NEW_ARRAY(
        __NLOLDVariable, __nlOLDCurrentContext->nb_variables
    );
    __nlOLDCurrentContext->alloc_variable = NLOLD_TRUE;
}

static void __nlOLDEndSystem() {
    __nlOLDTransition(__NLOLD_STATE_MATRIX_CONSTRUCTED, __NLOLD_STATE_SYSTEM_CONSTRUCTED);    
}

static void __nlOLDBeginMatrix() {
    NLOLDuint i;
    NLOLDuint n = 0;
    NLOLDenum storage = __NLOLD_ROWS;

    __nlOLDTransition(__NLOLD_STATE_SYSTEM, __NLOLD_STATE_MATRIX);

    for(i=0; i<__nlOLDCurrentContext->nb_variables; i++) {
        if(!__nlOLDCurrentContext->variable[i].locked) {
            __nlOLDCurrentContext->variable[i].index = n;
            n++;
        } else {
            __nlOLDCurrentContext->variable[i].index = ~0;
        }
    }

    __nlOLDCurrentContext->n = n;

    /* a least squares problem results in a symmetric matrix */
    if(__nlOLDCurrentContext->least_squares) {
        __nlOLDCurrentContext->symmetric = NLOLD_TRUE;
    }

    if(__nlOLDCurrentContext->symmetric) {
        storage = (storage | __NLOLD_SYMMETRIC);
    }

    /* SuperLU storage does not support symmetric storage */
    storage = (storage & ~__NLOLD_SYMMETRIC);

    __nlOLDSparseMatrixConstruct(&__nlOLDCurrentContext->M, n, n, storage);
    __nlOLDCurrentContext->alloc_M = NLOLD_TRUE;

    __nlOLDCurrentContext->x = __NLOLD_NEW_ARRAY(NLOLDfloat, n);
    __nlOLDCurrentContext->alloc_x = NLOLD_TRUE;
    
    __nlOLDCurrentContext->b = __NLOLD_NEW_ARRAY(NLOLDfloat, n);
    __nlOLDCurrentContext->alloc_b = NLOLD_TRUE;

    __nlOLDVariablesToVector();

    __nlOLDRowColumnConstruct(&__nlOLDCurrentContext->af);
    __nlOLDCurrentContext->alloc_af = NLOLD_TRUE;
    __nlOLDRowColumnConstruct(&__nlOLDCurrentContext->al);
    __nlOLDCurrentContext->alloc_al = NLOLD_TRUE;
    __nlOLDRowColumnConstruct(&__nlOLDCurrentContext->xl);
    __nlOLDCurrentContext->alloc_xl = NLOLD_TRUE;

    __nlOLDCurrentContext->current_row = 0;
}

static void __nlOLDEndMatrix() {
    __nlOLDTransition(__NLOLD_STATE_MATRIX, __NLOLD_STATE_MATRIX_CONSTRUCTED);    
    
    __nlOLDRowColumnDestroy(&__nlOLDCurrentContext->af);
    __nlOLDCurrentContext->alloc_af = NLOLD_FALSE;
    __nlOLDRowColumnDestroy(&__nlOLDCurrentContext->al);
    __nlOLDCurrentContext->alloc_al = NLOLD_FALSE;
    __nlOLDRowColumnDestroy(&__nlOLDCurrentContext->xl);
    __nlOLDCurrentContext->alloc_al = NLOLD_FALSE;
    
#if 0
    if(!__nlOLDCurrentContext->least_squares) {
        __nlOLD_assert(
            __nlOLDCurrentContext->current_row == 
            __nlOLDCurrentContext->n
        );
    }
#endif
}

static void __nlOLDBeginRow() {
    __nlOLDTransition(__NLOLD_STATE_MATRIX, __NLOLD_STATE_ROW);
    __nlOLDRowColumnZero(&__nlOLDCurrentContext->af);
    __nlOLDRowColumnZero(&__nlOLDCurrentContext->al);
    __nlOLDRowColumnZero(&__nlOLDCurrentContext->xl);
}

static void __nlOLDScaleRow(NLOLDfloat s) {
    __NLOLDRowColumn*    af = &__nlOLDCurrentContext->af;
    __NLOLDRowColumn*    al = &__nlOLDCurrentContext->al;
    NLOLDuint nf            = af->size;
    NLOLDuint nlOLD            = al->size;
    NLOLDuint i;
    for(i=0; i<nf; i++) {
        af->coeff[i].value *= s;
    }
    for(i=0; i<nlOLD; i++) {
        al->coeff[i].value *= s;
    }
    __nlOLDCurrentContext->right_hand_side *= s;
}

static void __nlOLDNormalizeRow(NLOLDfloat weight) {
    __NLOLDRowColumn*    af = &__nlOLDCurrentContext->af;
    __NLOLDRowColumn*    al = &__nlOLDCurrentContext->al;
    NLOLDuint nf            = af->size;
    NLOLDuint nlOLD            = al->size;
    NLOLDuint i;
    NLOLDfloat norm = 0.0;
    for(i=0; i<nf; i++) {
        norm += af->coeff[i].value * af->coeff[i].value;
    }
    for(i=0; i<nlOLD; i++) {
        norm += al->coeff[i].value * al->coeff[i].value;
    }
    norm = sqrt(norm);
    __nlOLDScaleRow(weight / norm);
}

static void __nlOLDEndRow() {
    __NLOLDRowColumn*    af = &__nlOLDCurrentContext->af;
    __NLOLDRowColumn*    al = &__nlOLDCurrentContext->al;
    __NLOLDRowColumn*    xl = &__nlOLDCurrentContext->xl;
    __NLOLDSparseMatrix* M  = &__nlOLDCurrentContext->M;
    NLOLDfloat* b        = __nlOLDCurrentContext->b;
    NLOLDuint nf          = af->size;
    NLOLDuint nlOLD          = al->size;
    NLOLDuint current_row = __nlOLDCurrentContext->current_row;
    NLOLDuint i;
    NLOLDuint j;
    NLOLDfloat S;
    __nlOLDTransition(__NLOLD_STATE_ROW, __NLOLD_STATE_MATRIX);

    if(__nlOLDCurrentContext->normalize_rows) {
        __nlOLDNormalizeRow(__nlOLDCurrentContext->row_scaling);
    } else {
        __nlOLDScaleRow(__nlOLDCurrentContext->row_scaling);
    }

    if(__nlOLDCurrentContext->least_squares) {
        for(i=0; i<nf; i++) {
            for(j=0; j<nf; j++) {
                __nlOLDSparseMatrixAdd(
                    M, af->coeff[i].index, af->coeff[j].index,
                    af->coeff[i].value * af->coeff[j].value
                );
            }
        }
        S = -__nlOLDCurrentContext->right_hand_side;
        for(j=0; j<nlOLD; j++) {
            S += al->coeff[j].value * xl->coeff[j].value;
        }
        for(i=0; i<nf; i++) {
            b[ af->coeff[i].index ] -= af->coeff[i].value * S;
        }
    } else {
        for(i=0; i<nf; i++) {
            __nlOLDSparseMatrixAdd(
                M, current_row, af->coeff[i].index, af->coeff[i].value
            );
        }
        b[current_row] = -__nlOLDCurrentContext->right_hand_side;
        for(i=0; i<nlOLD; i++) {
            b[current_row] -= al->coeff[i].value * xl->coeff[i].value;
        }
    }
    __nlOLDCurrentContext->current_row++;
    __nlOLDCurrentContext->right_hand_side = 0.0;    
    __nlOLDCurrentContext->row_scaling     = 1.0;
}

void nlOLDMatrixAdd(NLOLDuint row, NLOLDuint col, NLOLDfloat value)
{
    __NLOLDSparseMatrix* M  = &__nlOLDCurrentContext->M;
    __nlOLDCheckState(__NLOLD_STATE_MATRIX);
    __nlOLD_range_assert(row, 0, __nlOLDCurrentContext->n - 1);
    __nlOLD_range_assert(col, 0, __nlOLDCurrentContext->nb_variables - 1);
	__nlOLD_assert(!__nlOLDCurrentContext->least_squares);

	__nlOLDSparseMatrixAdd(M, row, col, value);
}

void nlOLDRightHandSideAdd(NLOLDuint index, NLOLDfloat value)
{
    NLOLDfloat* b = __nlOLDCurrentContext->b;

    __nlOLDCheckState(__NLOLD_STATE_MATRIX);
    __nlOLD_range_assert(index, 0, __nlOLDCurrentContext->n - 1);
	__nlOLD_assert(!__nlOLDCurrentContext->least_squares);

	b[index] += value;
}

void nlOLDCoefficient(NLOLDuint index, NLOLDfloat value) {
    __NLOLDVariable* v;
	unsigned int zero= 0;
    __nlOLDCheckState(__NLOLD_STATE_ROW);
    __nlOLD_range_assert(index, zero, __nlOLDCurrentContext->nb_variables - 1);
    v = &(__nlOLDCurrentContext->variable[index]);
    if(v->locked) {
        __nlOLDRowColumnAppend(&(__nlOLDCurrentContext->al), 0, value);
        __nlOLDRowColumnAppend(&(__nlOLDCurrentContext->xl), 0, v->value);
    } else {
        __nlOLDRowColumnAppend(&(__nlOLDCurrentContext->af), v->index, value);
    }
}

void nlOLDBegin(NLOLDenum prim) {
    switch(prim) {
    case NLOLD_SYSTEM: {
        __nlOLDBeginSystem();
    } break;
    case NLOLD_MATRIX: {
        __nlOLDBeginMatrix();
    } break;
    case NLOLD_ROW: {
        __nlOLDBeginRow();
    } break;
    default: {
        __nlOLD_assert_not_reached;
    }
    }
}

void nlOLDEnd(NLOLDenum prim) {
    switch(prim) {
    case NLOLD_SYSTEM: {
        __nlOLDEndSystem();
    } break;
    case NLOLD_MATRIX: {
        __nlOLDEndMatrix();
    } break;
    case NLOLD_ROW: {
        __nlOLDEndRow();
    } break;
    default: {
        __nlOLD_assert_not_reached;
    }
    }
}

/************************************************************************/
/* SuperLU wrapper */

/* Note: SuperLU is difficult to call, but it is worth it.    */
/* Here is a driver inspired by A. Sheffer's "cow flattener". */
static NLOLDboolean __nlOLDSolve_SUPERLU( NLOLDboolean do_perm) {

    /* OpenNLOLD Context */
    __NLOLDSparseMatrix* M  = &(__nlOLDCurrentContext->M);
    NLOLDfloat* b          = __nlOLDCurrentContext->b;
    NLOLDfloat* x          = __nlOLDCurrentContext->x;

    /* Compressed Row Storage matrix representation */
    NLOLDuint    n      = __nlOLDCurrentContext->n;
    NLOLDuint    nnz    = __nlOLDSparseMatrixNNZ(M); /* Number of Non-Zero coeffs */
    NLOLDint*    xa     = __NLOLD_NEW_ARRAY(NLOLDint, n+1);
    NLOLDfloat* rhs    = __NLOLD_NEW_ARRAY(NLOLDfloat, n);
    NLOLDfloat* a      = __NLOLD_NEW_ARRAY(NLOLDfloat, nnz);
    NLOLDint*    asub   = __NLOLD_NEW_ARRAY(NLOLDint, nnz);

    /* Permutation vector */
    NLOLDint*    perm_r  = __NLOLD_NEW_ARRAY(NLOLDint, n);
    NLOLDint*    perm    = __NLOLD_NEW_ARRAY(NLOLDint, n);

    /* SuperLU variables */
    SuperMatrix A, B; /* System       */
    SuperMatrix L, U; /* Inverse of A */
    NLOLDint info;       /* status code  */
    DNformat *vals = NULL; /* access to result */
    float *rvals  = NULL; /* access to result */

    /* SuperLU options and stats */
    superlu_options_t options;
    SuperLUStat_t     stat;


    /* Temporary variables */
    __NLOLDRowColumn* Ri = NULL;
    NLOLDuint         i,jj,count;
    
    __nlOLD_assert(!(M->storage & __NLOLD_SYMMETRIC));
    __nlOLD_assert(M->storage & __NLOLD_ROWS);
    __nlOLD_assert(M->m == M->n);
    
    
    /*
     * Step 1: convert matrix M into SuperLU compressed column 
     *   representation.
     * -------------------------------------------------------
     */

    count = 0;
    for(i=0; i<n; i++) {
        Ri = &(M->row[i]);
        xa[i] = count;
        for(jj=0; jj<Ri->size; jj++) {
            a[count]    = Ri->coeff[jj].value;
            asub[count] = Ri->coeff[jj].index;
            count++;
        }
    }
    xa[n] = nnz;

    /* Save memory for SuperLU */
    __nlOLDSparseMatrixClear(M);


    /*
     * Rem: symmetric storage does not seem to work with
     * SuperLU ... (->deactivated in main SLS::Solver driver)
     */
    sCreate_CompCol_Matrix(
        &A, n, n, nnz, a, asub, xa, 
        SLU_NR,              /* Row_wise, no supernode */
        SLU_S,               /* floats                */ 
        SLU_GE               /* general storage        */
    );

    /* Step 2: create vector */
    sCreate_Dense_Matrix(
        &B, n, 1, b, n, 
        SLU_DN, /* Fortran-type column-wise storage */
        SLU_S,  /* floats                          */
        SLU_GE  /* general                          */
    );
            

    /* Step 3: get permutation matrix 
     * ------------------------------
     * com_perm: 0 -> no re-ordering
     *           1 -> re-ordering for A^t.A
     *           2 -> re-ordering for A^t+A
     *           3 -> approximate minimum degree ordering
     */
    get_perm_c(do_perm ? 3 : 0, &A, perm);

    /* Step 4: call SuperLU main routine
     * ---------------------------------
     */

    set_default_options(&options);
    options.ColPerm = MY_PERMC;
    StatInit(&stat);

    sgssv(&options, &A, perm, perm_r, &L, &U, &B, &stat, &info);

    /* Step 5: get the solution
     * ------------------------
     * Fortran-type column-wise storage
     */
    vals = (DNformat*)B.Store;
    rvals = (float*)(vals->nzval);
    if(info == 0) {
        for(i = 0; i <  n; i++){
            x[i] = rvals[i];
        }
    }

    /* Step 6: cleanup
     * ---------------
     */

    /*
     *  For these two ones, onlOLDy the "store" structure
     * needs to be deallocated (the arrays have been allocated
     * by us).
     */
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);

    
    /*
     *   These ones need to be fully deallocated (they have been
     * allocated by SuperLU).
     */
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    StatFree(&stat);

    __NLOLD_DELETE_ARRAY(xa);
    __NLOLD_DELETE_ARRAY(rhs);
    __NLOLD_DELETE_ARRAY(a);
    __NLOLD_DELETE_ARRAY(asub);
    __NLOLD_DELETE_ARRAY(perm_r);
    __NLOLD_DELETE_ARRAY(perm);

    return (info == 0);
}


/************************************************************************/
/* nlOLDSolve() driver routine */

NLOLDboolean nlOLDSolve(void) {
    NLOLDboolean result = NLOLD_TRUE;

    __nlOLDCheckState(__NLOLD_STATE_SYSTEM_CONSTRUCTED);
    result = __nlOLDSolve_SUPERLU(NLOLD_TRUE);

    __nlOLDVectorToVariables();
    __nlOLDTransition(__NLOLD_STATE_SYSTEM_CONSTRUCTED, __NLOLD_STATE_SOLVED);

    return result;
}




