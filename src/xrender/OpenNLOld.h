/*
 *  $Id: ONL_opennl.h,v 1.3 2005/10/30 18:38:35 blendix Exp $
 *
 *  OpenNL: Numerical Library
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

/*
#define NL_DEBUG
#define NL_PARANOID
*/

#define NLOLD_USE_SUPERLU

#ifndef nloldOPENNL_H
#define nloldOPENNL_H

#ifdef __cplusplus
extern "C" {
#endif

#define NLOLD_VERSION_0_0 1

/* Datatypes */

typedef unsigned int	NLOLDenum;
typedef unsigned char	NLOLDboolean;
typedef unsigned int	NLOLDbitfield;
typedef void			NLOLDvoid;
typedef signed char		NLOLDbyte;		/* 1-byte signed */
typedef short			NLOLDshort;	/* 2-byte signed */
typedef int				NLOLDint;		/* 4-byte signed */
typedef unsigned char	NLOLDubyte;	/* 1-byte unsigned */
typedef unsigned short	NLOLDushort;	/* 2-byte unsigned */
typedef unsigned int	NLOLDuint;		/* 4-byte unsigned */
typedef int				NLOLDsizei;	/* 4-byte signed */
typedef float			NLOLDfloat;	/* single precision float */
typedef double			NLOLDdouble;	/* double precision float */

typedef void* NLOLDContext;

/* Constants */

#define NLOLD_FALSE   0x0
#define NLOLD_TRUE    0x1

/* Primitives */

#define NLOLD_SYSTEM  0x0
#define NLOLD_MATRIX  0x1
#define NLOLD_ROW     0x2

/* Solver Parameters */

#define NLOLD_SOLVER           0x100
#define NLOLD_NB_VARIABLES     0x101
#define NLOLD_LEAST_SQUARES    0x102
#define NLOLD_SYMMETRIC        0x106
#define NLOLD_ERROR            0x108

/* Enable / Disable */

#define NLOLD_NORMALIZE_ROWS  0x400

/* Row parameters */

#define NLOLD_RIGHT_HAND_SIDE 0x500
#define NLOLD_ROW_SCALING     0x501

/* Contexts */

NLOLDContext nlOLDNewContext(void);
void nlOLDDeleteContext(NLOLDContext context);
void nlOLDMakeCurrent(NLOLDContext context);
NLOLDContext nlOLDGetCurrent(void);

/* State get/set */

void nlOLDSolverParameterf(NLOLDenum pname, NLOLDfloat param);
void nlOLDSolverParameteri(NLOLDenum pname, NLOLDint param);

void nlOLDRowParameterf(NLOLDenum pname, NLOLDfloat param);
void nlOLDRowParameteri(NLOLDenum pname, NLOLDint param);

void nlOLDGetBooleanv(NLOLDenum pname, NLOLDboolean* params);
void nlOLDGetFloatv(NLOLDenum pname, NLOLDfloat* params);
void nlOLDGetIntergerv(NLOLDenum pname, NLOLDint* params);

void nlOLDEnable(NLOLDenum pname);
void nlOLDDisable(NLOLDenum pname);
NLOLDboolean nlOLDIsEnabled(NLOLDenum pname);

/* Variables */

void nlOLDSetVariable(NLOLDuint index, NLOLDfloat value);
NLOLDfloat nlOLDGetVariable(NLOLDuint index);
void nlOLDLockVariable(NLOLDuint index);
void nlOLDUnlockVariable(NLOLDuint index);
NLOLDboolean nlVariableIsLocked(NLOLDuint index);

/* Begin/End */

void nlOLDBegin(NLOLDenum primitive);
void nlOLDEnd(NLOLDenum primitive);
void nlOLDCoefficient(NLOLDuint index, NLOLDfloat value);

/* Setting random elements matrix/vector - not for least squares! */

void nlOLDMatrixAdd(NLOLDuint row, NLOLDuint col, NLOLDfloat value);
void nlOLDRightHandSideAdd(NLOLDuint index, NLOLDfloat value);

/* Solve */

NLOLDboolean nlOLDSolve(void);

#ifdef __cplusplus
}
#endif

#endif

