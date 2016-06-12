#include "stdafx.h"
#include "PackUV.h"
#include "BLI_edgehash.h"
#include "BLI_heap.h"
#include "MEM_guardedalloc.h"
#include "parametrizer_intern.h"
#include "parametrizer.h"


//static void hash_add_face(EdgeHash *ehash, MFace *mf)
//{
//	BLI_edgehash_insert(ehash, mf->v1, mf->v2, NULL);
//	BLI_edgehash_insert(ehash, mf->v2, mf->v3, NULL);
//	if(mf->v4) {
//		BLI_edgehash_insert(ehash, mf->v3, mf->v4, NULL);
//		BLI_edgehash_insert(ehash, mf->v4, mf->v1, NULL);
//	}
//	else
//		BLI_edgehash_insert(ehash, mf->v3, mf->v1, NULL);
//}

//void select_linked_tfaces_with_seams(int mode, unsigned int index,Mesh *gMeshs)			// Mesh *me,
//{
//	TFace *tf;
//	MFace *mf;
//	int a, doit=1, mark=0;
//	char *linkflag;
//	EdgeHash *ehash, *seamhash;
//	MEdge *med;
//	Mesh* me = gMeshs;
//
//	ehash= BLI_edgehash_new();
//	seamhash = BLI_edgehash_new();
//	linkflag= (char*)MEM_callocN(sizeof(char)*me->totface, "linkflaguv");
//
//	for(med=me->medge, a=0; a < me->totedge; a++, med++)
//		if(med->flag & ME_SEAM)
//			BLI_edgehash_insert(seamhash, med->v1, med->v2, NULL);
//
//	if (mode==0 || mode==1) {
//		/* only put face under cursor in array */
//		mf= ((MFace*)me->mface) + index;
//		hash_add_face(ehash, mf);
//		linkflag[index]= 1;
//	}
//	else {
//		/* fill array by selection */
//		tf= me->tface;
//		mf= me->mface;
//		for(a=0; a<me->totface; a++, tf++, mf++) {
//			if(tf->flag & TF_HIDE);
//			else if(tf->flag & TF_SELECT) {
//				hash_add_face(ehash, mf);
//				linkflag[a]= 1;
//			}
//		}
//	}
//
//	while(doit) {
//		doit= 0;
//
//		/* expand selection */
//		tf= me->tface;
//		mf= me->mface;
//		for(a=0; a<me->totface; a++, tf++, mf++) {
//			if(tf->flag & TF_HIDE)
//				continue;
//
//			if(!linkflag[a]) {
//				mark= 0;
//
//				if(!BLI_edgehash_haskey(seamhash, mf->v1, mf->v2))
//					if(BLI_edgehash_haskey(ehash, mf->v1, mf->v2))
//						mark= 1;
//				if(!BLI_edgehash_haskey(seamhash, mf->v2, mf->v3))
//					if(BLI_edgehash_haskey(ehash, mf->v2, mf->v3))
//						mark= 1;
//				if(mf->v4) {
//					if(!BLI_edgehash_haskey(seamhash, mf->v3, mf->v4))
//						if(BLI_edgehash_haskey(ehash, mf->v3, mf->v4))
//							mark= 1;
//					if(!BLI_edgehash_haskey(seamhash, mf->v4, mf->v1))
//						if(BLI_edgehash_haskey(ehash, mf->v4, mf->v1))
//							mark= 1;
//				}
//				else if(!BLI_edgehash_haskey(seamhash, mf->v3, mf->v1))
//					if(BLI_edgehash_haskey(ehash, mf->v3, mf->v1))
//						mark = 1;
//
//				if(mark) {
//					linkflag[a]= 1;
//					hash_add_face(ehash, mf);
//					doit= 1;
//				}
//			}
//		}
//
//	}
//
//	BLI_edgehash_free(ehash, NULL);
//	BLI_edgehash_free(seamhash, NULL);
//
//	if(mode==0 || mode==2) {
//		for(a=0, tf=me->tface; a<me->totface; a++, tf++)
//			if(linkflag[a])
//				tf->flag |= TF_SELECT;
//			else
//				tf->flag &= ~TF_SELECT;
//	}
//	else if(mode==1) {
//		for(a=0, tf=me->tface; a<me->totface; a++, tf++)
//			if(linkflag[a] && (tf->flag & TF_SELECT))
//				break;
//
//		if (a<me->totface) {
//			for(a=0, tf=me->tface; a<me->totface; a++, tf++)
//				if(linkflag[a])
//					tf->flag &= ~TF_SELECT;
//		}
//		else {
//			for(a=0, tf=me->tface; a<me->totface; a++, tf++)
//				if(linkflag[a])
//					tf->flag |= TF_SELECT;
//		}
//	}
//
//	MEM_freeN(linkflag);
//
//	//	BIF_undo_push("Select linked UV face");
//	//	object_tface_flags_changed(OBACT, 0);
//}

ParamHandle *construct_new_param_handle(ParamHandle * old_handle, vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<int> &visited, bool fill, bool implicit)
{

	int  m_NumberOfFaces = pFaces.size();
	int  m_NumberOfVertices = pVertices.size();

	ParamHandle *handle;
	handle = param_construct_begin();

	PChart *chart;
	PHandle *phandle = (PHandle*)old_handle;

	for (int i = 0; i < phandle->ncharts; i++) 
		{
			chart = phandle->charts[i];
			PFace *f;
			for (f=chart->faces; f; f=f->nextlink) {


				ParamKey key, vkeys[4];
				bool pin[4], select[4];							// was ParamBool
				float *co[4];
				float *uv[4];
				int nverts;

				int a = f->idx;

				key = (ParamKey)a;
				vkeys[0] = (ParamKey)pFaces[a].v[0];
				vkeys[1] = (ParamKey)pFaces[a].v[1];
				vkeys[2] = (ParamKey)pFaces[a].v[2];

				co[0] = pVertices[pFaces[a].v[0]].getAddress();
				co[1] = pVertices[pFaces[a].v[1]].getAddress();
				co[2] = pVertices[pFaces[a].v[2]].getAddress();

				uv[0] = NULL;
				uv[1] = NULL;
				uv[2] = NULL;

				pin[0] = false;
				pin[1] = false;
				pin[2] = false;

				select[0] = false;
				select[1] = false;
				select[2] = false;

				//if (mf->v4) {
				//	vkeys[3] = (ParamKey)mf->v4;
				//	co[3] = (mv+mf->v4)->co;
				//	uv[3] = tf->uv[3];
				//	pin[3] = false;
				//	select[3] = ((tf->flag & TF_SEL4) != 0);
				//	nverts = 4;
				//}
				//else
				nverts = 3;

				param_face_add(handle, key, nverts, vkeys, co, uv, pin, select);


			}
	}

	for (int i = 0; i < phandle->ncharts; i++) 
		{
			chart = phandle->charts[i];
			PEdge *e;

			for (e = chart->edges; e; e = e->nextlink) {

				PFace* f = e->face;
				if (e->pair) {
					PFace* pair_f = e->pair->face;

					int f_index = e->face->idx;
					int pair_f_index = e->pair->face->idx;

					if (visited[f_index] != visited[pair_f_index]) {

						ParamKey vkeys[2];
						vkeys[0] = (ParamKey)e->vert->u.key;
						vkeys[1] = (ParamKey)e->pair->vert->u.key;
						param_edge_set_seam(handle, vkeys);

						
					}

				}

			}
	}



	param_construct_end(handle, fill != 0, implicit != 0);

	return handle;

}

ParamHandle *construct_param_handle_from_data(vector <VEC3>  &pVertices, vector<CFace> &pFaces, bool fill, bool implicit)
{
	int  m_NumberOfFaces = pFaces.size();
	int  m_NumberOfVertices = pVertices.size();

	ParamHandle *handle;
	handle = param_construct_begin();

	for (int a=0; a<m_NumberOfFaces; a++) {
		ParamKey key, vkeys[4];
		bool pin[4], select[4];							// was ParamBool
		float *co[4];
		float *uv[4];
		int nverts;


		key = (ParamKey)a;
		vkeys[0] = (ParamKey)pFaces[a].v[0];
		vkeys[1] = (ParamKey)pFaces[a].v[1];
		vkeys[2] = (ParamKey)pFaces[a].v[2];

		co[0] = pVertices[pFaces[a].v[0]].getAddress();
		co[1] = pVertices[pFaces[a].v[1]].getAddress();
		co[2] = pVertices[pFaces[a].v[2]].getAddress();

		uv[0] = NULL;
		uv[1] = NULL;
		uv[2] = NULL;

		pin[0] = false;
		pin[1] = false;
		pin[2] = false;

		select[0] = false;
		select[1] = false;
		select[2] = false;

		//if (mf->v4) {
		//	vkeys[3] = (ParamKey)mf->v4;
		//	co[3] = (mv+mf->v4)->co;
		//	uv[3] = tf->uv[3];
		//	pin[3] = false;
		//	select[3] = ((tf->flag & TF_SEL4) != 0);
		//	nverts = 4;
		//}
		//else
		nverts = 3;

		param_face_add(handle, key, nverts, vkeys, co, uv, pin, select);
	}


	param_construct_data_end(handle, fill != 0, implicit != 0);

	return handle;

}


ParamHandle *construct_param_handle(vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<int> &visited, bool fill, bool implicit)
{

	int  m_NumberOfFaces = pFaces.size();
	int  m_NumberOfVertices = pVertices.size();

	ParamHandle *handle;
	handle = param_construct_begin();

	for (int a=0; a<m_NumberOfFaces; a++) {
		ParamKey key, vkeys[4];
		bool pin[4], select[4];							// was ParamBool
		float *co[4];
		float *uv[4];
		int nverts;


		key = (ParamKey)a;
		vkeys[0] = (ParamKey)pFaces[a].v[0];
		vkeys[1] = (ParamKey)pFaces[a].v[1];
		vkeys[2] = (ParamKey)pFaces[a].v[2];

		co[0] = pVertices[pFaces[a].v[0]].getAddress();
		co[1] = pVertices[pFaces[a].v[1]].getAddress();
		co[2] = pVertices[pFaces[a].v[2]].getAddress();

		uv[0] = NULL;
		uv[1] = NULL;
		uv[2] = NULL;

		pin[0] = false;
		pin[1] = false;
		pin[2] = false;

		select[0] = false;
		select[1] = false;
		select[2] = false;

		//if (mf->v4) {
		//	vkeys[3] = (ParamKey)mf->v4;
		//	co[3] = (mv+mf->v4)->co;
		//	uv[3] = tf->uv[3];
		//	pin[3] = false;
		//	select[3] = ((tf->flag & TF_SEL4) != 0);
		//	nverts = 4;
		//}
		//else
		nverts = 3;

		param_face_add(handle, key, nverts, vkeys, co, uv, pin, select);
	}

	PHandle *phandle = (PHandle*)handle;
	phandle->chart_index = &visited;


		
	param_construct_data_end(handle, fill != 0, implicit != 0);

	return handle;

}







//
