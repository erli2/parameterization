// sample.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <time.h>

int _tmain(int argc, _TCHAR* argv[])
{


		char* OutFileName = "..\\data\\desk.obj";
		FILE* OutFileHandle = fopen(OutFileName,"r");

		vector<VEC3>vertex;
	    vector<CFace>face;

		const int MAX_LINE= 100;
		char strBuff[MAX_LINE]	= {0};
		char chKeyword          = 0;

		while(EOF != fscanf_s(OutFileHandle, "%s", strBuff, MAX_LINE)){
			chKeyword = strBuff[0];

			if(chKeyword == 'v' && strBuff[1]!='t'){

				VEC3 vNewVertex;
				fscanf_s(OutFileHandle, "%f %f %f", &vNewVertex.x, &vNewVertex.y, &vNewVertex.z);
				vertex.push_back(vNewVertex);

			}
			if(chKeyword == 'f'){

				CFace vNewFace;
				int x,y,z;
				fscanf_s(OutFileHandle, "%d/%d %d/%d %d/%d", &vNewFace.v[0],&x, &vNewFace.v[1], &y,&vNewFace.v[2],&z);
				vNewFace.v[0] -= 1;
				vNewFace.v[1] -= 1;
				vNewFace.v[2] -= 1;
				face.push_back(vNewFace);
			}
		}

		vector<VEC2> pUVs;
		vector<CFace> UVIndex;

		clock_t start,end;
		start = clock();
		
  		//UVfun::AutoMapping(vertex, face, pUVs, UVIndex,false,60.0);
		UVfun::BoxMapping(vertex, face, pUVs, UVIndex);

		end = clock();
		double dur = (double)(end - start);
		printf("Use Time:%f\n",(dur/CLOCKS_PER_SEC));

		return 0;

	
}

