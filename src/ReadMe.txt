========================================================================
     Auto mapping for 3d triangular meshes
========================================================================

xrender 
    This is for mesh segmentation and parameterization.
    
superLU 
    This is linear solver.

sample 
    This is for tesing the algorithms.
    

/////////////////////////////////////////////////////////////////////////////
Usage 
    Please see sample.cpp file for details.
    

	 ------------------------------------
	  AutoMapping function
	 ------------------------------------
	 --- pVertices： vertices of the input mesh 
	 --- pFaces：faces of the input mesh 
	 --- pUVs: output for the UV texture coordinate 
	 --- UVIndex：output for the UV indices of each face 
	 --- use_LSCM：using Least Square Comformal Mapping when it is True 
	               
	 --- angle：angle of auto flatten mapping 
	 --- margin: maringal threshold for the packing 

	static void AutoMapping(vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<VEC2> &pUVs,vector<CFace>& UVIndex, bool use_LSCM = false, float angle = 45.0,float margin = 0.02);

/////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////
