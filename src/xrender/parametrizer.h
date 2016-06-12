#ifndef __PARAMETRIZER_H__
#define __PARAMETRIZER_H__

#ifdef __cplusplus
extern "C" {
#endif


typedef void ParamHandle;	/* handle to a set of charts */
typedef long ParamKey;		/* (hash) key for identifying verts and faces */
typedef enum ParamBool {
	PARAM_TRUE = 1,
	PARAM_FALSE = 0
} ParamBool;

/* Chart construction:
   -------------------
   - faces and seams may only be added between construct_{begin|end}
   - the pointers to co and uv are stored, rather than being copied
   - vertices are implicitly created
   - in construct_end the mesh will be split up according to the seams
   - the resulting charts must be:
      - manifold, connected, open (at least one boundary loop)
   - output will be written to the uv pointers
*/

#define MAX2(x,y)               ( (x)>(y) ? (x) : (y) )
#define MAX3(x,y,z)             MAX2( MAX2((x),(y)) , (z) )
#define INIT_MINMAX2(min, max) { (min)[0]= (min)[1]= 1.0e30f; (max)[0]= (max)[1]= -1.0e30f; }
#define DO_MINMAX(vec, min, max) { if( (min)[0]>(vec)[0] ) (min)[0]= (vec)[0]; \
	if( (min)[1]>(vec)[1] ) (min)[1]= (vec)[1]; \
	if( (min)[2]>(vec)[2] ) (min)[2]= (vec)[2]; \
	if( (max)[0]<(vec)[0] ) (max)[0]= (vec)[0]; \
	if( (max)[1]<(vec)[1] ) (max)[1]= (vec)[1]; \
	if( (max)[2]<(vec)[2] ) (max)[2]= (vec)[2]; } \

#define DO_MINMAX2(vec, min, max) { if( (min)[0]>(vec)[0] ) (min)[0]= (vec)[0]; \
	if( (min)[1]>(vec)[1] ) (min)[1]= (vec)[1]; \
	if( (max)[0]<(vec)[0] ) (max)[0]= (vec)[0]; \
	if( (max)[1]<(vec)[1] ) (max)[1]= (vec)[1]; }

ParamHandle *param_construct_begin();

void param_face_add(ParamHandle *handle,
                    ParamKey key,
                    int nverts,	
                    ParamKey *vkeys,
                    float **co,
                    float **uv,
					bool *pin,
					bool *select);

void param_edge_set_seam(ParamHandle *handle,
                         ParamKey *vkeys);

void param_construct_data_end(ParamHandle *handle, bool fill, bool impl);

void param_construct_end(ParamHandle *handle, bool fill, bool impl);

void param_delete(ParamHandle *chart);

/* Least Squares Conformal Maps:
   -----------------------------
   - charts with less than two pinned vertices are assigned 2 pins
   - lscm is divided in three steps:
      - begin: compute matrix and it's factorization (expensive)
      - solve using pinned coordinates (cheap)
	  - end: clean up 
	- uv coordinates are allowed to change within begin/end, for
	  quick re-solving
*/

void param_lscm_begin(ParamHandle *handle, bool live, bool abf);
void param_lscm_solve(ParamHandle *handle);
void param_lscm_end(ParamHandle *handle);

/* Stretch */

void param_stretch_begin(ParamHandle *handle);
void param_stretch_blend(ParamHandle *handle, float blend);
void param_stretch_iter(ParamHandle *handle);
void param_stretch_end(ParamHandle *handle);

/* Area Smooth */

void param_smooth_area(ParamHandle *handle);

/* Packing */

//void param_pack(ParamHandle *handle);
void param_pack(ParamHandle *handle, float margin);

/* Flushing */

void param_flush(ParamHandle *handle);
void param_flush_restore(ParamHandle *handle);


#ifdef __cplusplus
}
#endif

#endif /*__PARAMETRIZER_H__*/
