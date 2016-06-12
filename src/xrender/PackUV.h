#ifndef __PACKUV_H__
#define __PACKUV_H__
#include "AutoFlattenMapping.h"
#include "parametrizer.h"
#include "BLI_memarena.h"

//void pack_UV(float margin, Mesh* gMeshes);
//void select_linked_tfaces_with_seams(int mode, unsigned int index,Mesh* gMeshs);
ParamHandle *construct_param_handle_from_data(vector <VEC3>  &pVertices, vector<CFace> &pFaces, bool fill = false, bool implicit = false);
ParamHandle *construct_new_param_handle(ParamHandle * old_handle, vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<int> &visited, bool fill = false, bool implicit = false);
ParamHandle *construct_param_handle(vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<int> &visited, bool fill = false, bool implicit = false);
#endif