/*
###############################################################
Created on Sat May 11 16:07:02 2015
:author: Er Li {lier198631@126.com}
###############################################################
*/

#ifndef __AUTOFLATTENMAPPING_H__
#define __AUTOFLATTENMAPPING_H__

#define DDF_DEBUG 1

#include <assert.h>
#include "declare.h"
#include <vector>
#include <map>
#include "vector.h"

using namespace std;

struct CFace
{
	int v[3];
};

class VEC2
{
public:
	union
	{
		struct 
		{
			float x, y;//!<用两个Real类型变量来表示一个二维行向量  
		};
		struct 
		{
			float u, v;//!<用两个Real类型变量来表示一个二维行向量
		};
		float m[2];//!<用一个数组来表示一个二维行向量
	};

	VEC2 (const float x, const float y)
	{
		this->x = x;
		this->y = y;
	}
};

typedef BasicVector::Vector3Dim<float> VEC3;



class  DEVICE_API UVfun
{
public:
	UVfun();
	~UVfun();

	/* 自动展UV函数 AutoMapping
	 * pVertices：输入顶点坐标； pFaces：输入面片信息； pUVs:输出UV坐标；UVIndex：输出mesh每个面对应的UV索引；angle：展平角度，默认 45度,最大90度，角度越小，分割越破碎； margin:紧缩阈值，范围0.0 至1.0,默认0.02；
	 */

	static void AutoMapping(vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<VEC2> &pUVs,vector<CFace>& UVIndex, bool use_LSCM = false, float angle = 45.0,float margin = 0.02);


	/* 自动展UV函数 BoxMapping
	 * pVertices：输入顶点坐标； pFaces：输入面片信息； pUVs:输出UV坐标；UVIndex：输出mesh每个面对应的UV索引；margin:紧缩阈值，范围0.0 至1.0,默认0.02；
	 */

	static void BoxMapping(vector <VEC3>  &pVertices, vector<CFace> &pFaces, vector<VEC2> &pUVs,vector<CFace>& UVIndex, bool use_LSCM = false, float margin = 0.02);
	


};



#endif
