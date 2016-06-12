/*
###############################################################
Created on Sat May 11 16:07:02 2015
:author: Er Li {lier198631@126.com}
###############################################################
*/

#include "stdafx.h"
#include "AutoFlattenMapping.h "
#include "PackUV.h"
#include "BLI_memarena.h"
#include "BLI_heap.h"
#include "parametrizer_intern.h"

#include<list>

// Perfrom planar projection for each chart
void ChartProjection(ParamHandle *handle, vector<VEC3> &chart_normals, vector<int> &visited)
{
		PChart *chart;
		PHandle *phandle = (PHandle*)handle;

		VEC3  positive_x_direction(1.0,0.0,0.0);
		VEC3  positive_y_direction(0.0,1.0,0.0);


		for (int i = 0; i < phandle->ncharts; i++) 
		{
			chart = phandle->charts[i];
			PFace *f;
			int index = chart->faces->idx;
			int chart_indexs = visited[index];	


			VEC3 Nchart = chart_normals[chart_indexs];
			VEC3 Nx,Ny;
			Nx = BasicVector::cross(positive_x_direction,Nchart);
			if(BasicVector::norm(Nx)<10e-5)
				Nx = positive_y_direction;
			Nx = BasicVector::getNormalized(Nx);
		
			Ny = BasicVector::cross(Nx,Nchart);
			Ny = BasicVector::getNormalized(Ny);


			for (f=chart->faces; f; f=f->nextlink) {

				int index = f->idx;
				int chart_index = visited[index];
				if (chart_index != chart_indexs)
					break;

				PEdge *e1 = f->edge, *e2 = e1->next, *e3 = e2->next;
				PVert *v1 = e1->vert, *v2 = e2->vert, *v3 = e3->vert;

				VEC3 v_vector =VEC3(v1->co) ;
				float Ux = BasicVector::dot(v_vector,Nx);
				float Uy = BasicVector::dot(v_vector,Ny);
				if (!(v1->flag & PVERT_VISIT)) {
					v1->uv[0] = Ux;
					v1->uv[1] = Uy;
					v1->flag |= PVERT_VISIT;
				}


				v_vector =VEC3(v2->co) ;
				Ux = BasicVector::dot(v_vector,Nx);
				Uy = BasicVector::dot(v_vector,Ny);
				if (!(v2->flag & PVERT_VISIT)) {
					v2->uv[0] = Ux;
					v2->uv[1] = Uy;
					v2->flag |= PVERT_VISIT;
				}

				v_vector =VEC3(v3->co) ;
				Ux = BasicVector::dot(v_vector,Nx);
				Uy = BasicVector::dot(v_vector,Ny);
				if (!(v3->flag & PVERT_VISIT)) {
					v3->uv[0] = Ux;
					v3->uv[1] = Uy;
					v3->flag |= PVERT_VISIT;
				}
			}
		}
}

//calculate normal for each face
void CalculateNormals(ParamHandle *handle, vector<VEC3> &normals)
{
		PChart *chart;
		PHandle *phandle = (PHandle*)handle;

		for (int i = 0; i < phandle->ncharts; i++) 
		{
			chart = phandle->charts[i];
			PFace *f;
			
			for (f=chart->faces; f; f=f->nextlink) {

				
				PEdge *e1 = f->edge, *e2 = e1->next, *e3 = e2->next;
				PVert *v1 = e1->vert, *v2 = e2->vert, *v3 = e3->vert;

				VEC3 v1v0 = VEC3(v2->co) - VEC3(v1->co);
				VEC3 v2v1 = VEC3(v3->co) - VEC3(v2->co);

				VEC3 normal = BasicVector::cross(v1v0,v2v1);
				normal = BasicVector::getNormalized(normal);

				int index = f->idx;
				normals[index] = normal;
				
			}
		}
}


void AssignBoxChart(vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<VEC3> &normals, vector<int>& visited, VEC3 *xyz_axis)
{
		for (int i = 0; i < pFaces.size(); i++) 
		{

				int vindex[3];

				vindex[0] = (ParamKey)pFaces[i].v[0];
				vindex[1] = (ParamKey)pFaces[i].v[1];
				vindex[2] = (ParamKey)pFaces[i].v[2];

				VEC3 v1v0 = pVertices[vindex[1]] - pVertices[vindex[0]];
				VEC3 v2v1 = pVertices[vindex[2]] - pVertices[vindex[1]];

				VEC3 normal = BasicVector::cross(v1v0,v2v1);
				normal = BasicVector::getNormalized(normal);
				
				normals[i] = normal;
				double max_variation = -100.0;
				int chart_index = 0;

				for (int j = 0; j < 6; j++){
					double angle_variation = BasicVector::dot(normal, xyz_axis[j]);
					if (angle_variation>max_variation){
						max_variation = angle_variation;
						chart_index = j;
					}
				}
				visited[i] = chart_index + 1;
				
		}
}


/* Get the face from which we grow the chart */
PFace* GetStartFace(ParamHandle *handle, vector<VEC3> &normals, vector<int>& visited,  int m_NumberOfFaces, VEC3 normal,float thresh)
{

		PChart *chart;
		int i;
		PHandle *phandle = (PHandle*)handle;

		PFace* s = NULL;
		double max = -100.0;

		for (i = 0; i < phandle->ncharts; i++) 
		{
			chart = phandle->charts[i];
			PFace *f;

			for (f=chart->faces; f; f=f->nextlink) {

				int index = f->idx;
				if (visited[index])
					continue;
				double t = BasicVector::dot(normal,normals[index]);
				if ( t > thresh && t > max) {
					s = f;
					max = t;
				}

			}
		}

		return s;

}


	
void UVfun::AutoMapping(vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<VEC2> &pUVs,vector<CFace>& UVIndex, bool use_LSCM, float angle, float margin)
	{
		ParamHandle *handle = construct_param_handle_from_data(pVertices, pFaces, false,  0);

		int  m_NumberOfFaces = pFaces.size();
		vector<int>visited(m_NumberOfFaces,0);

		vector<VEC3>normals(m_NumberOfFaces);

		UVIndex.resize(m_NumberOfFaces);

		CalculateNormals(handle, normals);

		VEC3  x_axis(1.0,0.0,0.0);
		VEC3  y_axis(0.0,1.0,0.0);
		VEC3  z_axis(0.0,0.0,1.0);

		VEC3 xyz_axis[6];
		xyz_axis[0] = x_axis;
		xyz_axis[1] = x_axis*-1.0;
		xyz_axis[2] = y_axis;
		xyz_axis[3] = y_axis*-1.0;
		xyz_axis[4] = z_axis;
		xyz_axis[5] = z_axis*-1.0;

		int d = 0;
		VEC3 normal;
		PFace* start_face = NULL;

		float thresh = float(cos(2* M_PI*angle/360.0));
				
		while (d < 6) {
			start_face = GetStartFace(handle, normals, visited, m_NumberOfFaces, xyz_axis[d], thresh);
			if (start_face != NULL) {
				normal = xyz_axis[d];
				d++;
				break;
			}
			d++;
		}

		if (start_face == NULL) {
			start_face = GetStartFace(handle, normals, visited, m_NumberOfFaces, xyz_axis[0], -10.0);
			normal = normals[start_face->idx];
		}

		std::vector<VEC3> chart_normal;//(chart_num,VEC3(0.0,0.0,0.0));
		chart_normal.push_back(VEC3(0.0,0.0,0.0));
		chart_normal.push_back(normal);

		int chart_num = 1;
		
		while (start_face != NULL) {

			int index = start_face->idx;
			if (visited[index])
				continue;
			visited[index] = chart_num;
			normal = BasicVector::getNormalized(normal);

			std::list<PFace*> faces;
			faces.push_back(start_face);

			while (!faces.empty()) {
				PFace *f = faces.front();
				faces.pop_front();
				int index_f = f->idx;

				PEdge* e = f->edge;
				do {

					PEdge* pair_e = e->pair;
					e = e->next;
					PEdge* we = pair_e;
					
					if (pair_e) {
						
						PFace* pair_f = pair_e->face;
						int pair_f_index = pair_f->idx;
						if (visited[pair_f_index]) 
							continue;

						VEC3 normal_f =normals[pair_f_index];

						if (BasicVector::dot(normal_f, normal) > thresh) {

							faces.push_front(pair_f);
							visited[pair_f_index] = chart_num;
						}

					}
					

					
				} while ( e && e != f->edge);

			}

			if ( d == 6 ) d = 0;

			PFace* new_start_face = NULL;

			while (d < 6) {
				new_start_face = GetStartFace(handle, normals, visited, m_NumberOfFaces, xyz_axis[d], thresh);
				if (new_start_face != NULL) {
					normal = xyz_axis[d];
					d++;
					break;
				}
				d++;
		    }

			if (d == 6 && new_start_face == NULL) {
				new_start_face = GetStartFace(handle, normals, visited, m_NumberOfFaces, xyz_axis[0], -10.0);
				if (new_start_face != NULL) {
					normal = normals[new_start_face->idx];
				}
			}

				start_face = new_start_face;
				chart_normal.push_back(normal);
				chart_num++;
		}


		std::vector<std::vector<int> > chart_index;
		chart_index.resize(chart_num);

		for(int index = 0; index < m_NumberOfFaces; index++ ) {
			chart_index[visited[index]].push_back(index);
		}	

		for(int index = 1; index <chart_num; index++ ) {
			chart_normal[index] = BasicVector::getNormalized(chart_normal[index]);		
		}



		ParamHandle *new_handle = construct_new_param_handle(handle,pVertices, pFaces, visited);
		
		if(!use_LSCM)
			ChartProjection(new_handle, chart_normal, visited);

		else{

			param_lscm_begin(new_handle, PARAM_FALSE, true);
			param_lscm_solve(new_handle);
			param_lscm_end(new_handle);
		}
		//param_smooth_area(handle);
		param_pack(new_handle,margin);

		param_flush(new_handle);


		/* This is for test */

		char* OutFileName = "..\\data\\test.obj";
		FILE* OutFileHandle = fopen(OutFileName,"w");

		
		int index = 0;

		PChart *chart;
		PHandle *phandle = (PHandle*)new_handle;


		for (int i = 0; i < phandle->ncharts; i++) {
			chart = phandle->charts[i];

			PVert *v;

			for (v=chart->verts; v; v=v->nextlink) {

				if (v->uv) {
					VEC2 tc(v->uv[0],v->uv[1]);

					pUVs.push_back(tc);
					v->u.id = index;
					index++;

				}
			}


		}

		// for test
		for(unsigned int i = 0; i < pVertices.size(); i++)
		{
			float x = pVertices[i].x;
			float y = pVertices[i].y;
			float z = pVertices[i].z;
			fprintf( OutFileHandle, "v %f %f %f\n", x, y, z); //use for test
		}

		for(unsigned int i = 0; i < pUVs.size(); i++)
		{
			float x = pUVs[i].x;
			float y = pUVs[i].y;
			fprintf( OutFileHandle, "vt %f %f\n", x, y);
		}


		for (int i = 0; i < phandle->ncharts; i++) {
			chart = phandle->charts[i];

			PFace *f;

			for (f=chart->faces; f; f=f->nextlink) {

					PEdge *e1 = f->edge, *e2 = e1->next, *e3 = e2->next;
					PVert *v1 = e1->vert, *v2 = e2->vert, *v3 = e3->vert;

					int V1Index = (v1->co - pVertices[0].getAddress())/3 + 1;
					int V2Index = (v2->co - pVertices[0].getAddress())/3 + 1;
					int V3Index = (v3->co - pVertices[0].getAddress())/3 + 1;


					int UV1Index = v1->u.id + 1 ;
					int UV2Index = v2->u.id + 1 ;
					int UV3Index = v3->u.id + 1 ;

					int a = f->idx;

					//int V1Index = pFaces[a].v[0] + 1;
					//int V2Index = pFaces[a].v[1] + 1;
					//int V3Index = pFaces[a].v[2] + 1;

					CFace vNewFace;

					vNewFace.v[0] = UV1Index - 1;
					vNewFace.v[1] = UV2Index - 1;
					vNewFace.v[2] = UV3Index - 1;

					UVIndex[a] = vNewFace;

					fprintf(OutFileHandle, "f %d/%d %d/%d %d/%d\n",V1Index,UV1Index,V2Index,UV2Index,V3Index,UV3Index); //use for test
					
			}
			}


				
		fclose ( OutFileHandle );

		param_delete(handle);
		param_delete(new_handle);
		

	}


	void UVfun::BoxMapping(vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<VEC2> &pUVs,vector<CFace>& UVIndex, bool use_LSCM, float margin)
	{
		
		int  m_NumberOfFaces = pFaces.size();
		vector<int>visited(m_NumberOfFaces,0);

		vector<VEC3>normals(m_NumberOfFaces);

		UVIndex.resize(m_NumberOfFaces);
				

		VEC3  x_axis(1.0,0.0,0.0);
		VEC3  y_axis(0.0,1.0,0.0);
		VEC3  z_axis(0.0,0.0,1.0);

		VEC3 xyz_axis[6];
		xyz_axis[0] = x_axis;
		xyz_axis[1] = x_axis*-1.0;
		xyz_axis[2] = y_axis;
		xyz_axis[3] = y_axis*-1.0;
		xyz_axis[4] = z_axis;
		xyz_axis[5] = z_axis*-1.0;

		int chart_num = 7;

		std::vector<VEC3> chart_normal(chart_num);

		AssignBoxChart(pVertices, pFaces, normals, visited, xyz_axis);

		for(int index = 1; index <chart_num; index++ ) {
			chart_normal[index] = xyz_axis[index-1];		
		}

		
		ParamHandle *handle = construct_param_handle(pVertices, pFaces, visited);
		
		if(!use_LSCM)
			ChartProjection(handle, chart_normal, visited);

		else{

			param_lscm_begin(handle, PARAM_FALSE, true);
			param_lscm_solve(handle);
			param_lscm_end(handle);
		}
		//param_smooth_area(handle);
		param_pack(handle,margin);

		param_flush(handle);


		/* This is for test */

		char* OutFileName = "..\\data\\test.obj";
		FILE* OutFileHandle = fopen(OutFileName,"w");

		
		int index = 0;

		PChart *chart;
		PHandle *phandle = (PHandle*)handle;


		for (int i = 0; i < phandle->ncharts; i++) {
			chart = phandle->charts[i];

			PVert *v;

			for (v=chart->verts; v; v=v->nextlink) {

				if (v->uv) {
					VEC2 tc(v->uv[0],v->uv[1]);

					pUVs.push_back(tc);
					v->u.id = index;
					index++;

				}
			}


		}

		// for test
		for(unsigned int i = 0; i < pVertices.size(); i++)
		{
			float x = pVertices[i].x;
			float y = pVertices[i].y;
			float z = pVertices[i].z;
			fprintf( OutFileHandle, "v %f %f %f\n", x, y, z); //use for test
		}

		for(unsigned int i = 0; i < pUVs.size(); i++)
		{
			float x = pUVs[i].x;
			float y = pUVs[i].y;
			fprintf( OutFileHandle, "vt %f %f\n", x, y);
		}


		for (int i = 0; i < phandle->ncharts; i++) {
			chart = phandle->charts[i];

			PFace *f;

			for (f=chart->faces; f; f=f->nextlink) {

					PEdge *e1 = f->edge, *e2 = e1->next, *e3 = e2->next;
					PVert *v1 = e1->vert, *v2 = e2->vert, *v3 = e3->vert;

					int V1Index = (v1->co - pVertices[0].getAddress())/3 + 1;
					int V2Index = (v2->co - pVertices[0].getAddress())/3 + 1;
					int V3Index = (v3->co - pVertices[0].getAddress())/3 + 1;


					int UV1Index = v1->u.id + 1 ;
					int UV2Index = v2->u.id + 1 ;
					int UV3Index = v3->u.id + 1 ;

					int a = f->idx;

					//int V1Index = pFaces[a].v[0] + 1;
					//int V2Index = pFaces[a].v[1] + 1;
					//int V3Index = pFaces[a].v[2] + 1;

					CFace vNewFace;

					vNewFace.v[0] = UV1Index - 1;
					vNewFace.v[1] = UV2Index - 1;
					vNewFace.v[2] = UV3Index - 1;

					UVIndex[a] = vNewFace;

					fprintf(OutFileHandle, "f %d/%d %d/%d %d/%d\n",V1Index,UV1Index,V2Index,UV2Index,V3Index,UV3Index); //use for test
					
			}
			}


				
		fclose ( OutFileHandle );

		param_delete(handle);
		//param_delete(new_handle);
		

	}




	

