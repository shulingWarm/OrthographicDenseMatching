////////////////////////////////////////////////////////////////////////////
//	File:		SiftMatch.cpp
//	Author:		Changchang Wu
//	Description :	implementation of SiftMatchGPU and SiftMatchGL
//
//
//	Copyright (c) 2007 University of North Carolina at Chapel Hill
//	All Rights Reserved
//
//	Permission to use, copy, modify and distribute this software and its
//	documentation for educational, research and non-profit purposes, without
//	fee, and without a written agreement is hereby granted, provided that the
//	above copyright notice and the following paragraph appear in all copies.
//	
//	The University of North Carolina at Chapel Hill make no representations
//	about the suitability of this software for any purpose. It is provided
//	'as is' without express or implied warranty. 
//
//	Please send BUG REPORTS to ccwu@cs.unc.edu
//
////////////////////////////////////////////////////////////////////////////


#include "GL/glew.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <strstream>
#include <algorithm>
using namespace std;
#include <string.h>
#include "GlobalUtil.h"

#if defined(CG_SIFTGPU_ENABLED) 
#include "ProgramCG.h"
#endif

#include "ProgramGLSL.h"
#include "GLTexImage.h"
#include "SiftGPU.h"
#include "SiftMatch.h"
#include "FrameBufferObject.h"

#if defined(CUDA_SIFTGPU_ENABLED)
#include "CuTexImage.h"
#include "SiftMatchCU.h"
#endif


SiftMatchGL::SiftMatchGL(int max_sift, int use_glsl): SiftMatchGPU()
{
	s_multiply = s_col_max = s_row_max = s_guided_mult = NULL;
	_num_sift[0] = _num_sift[1] = 0;
	_id_sift[0] = _id_sift[1] = 0;
	_have_loc[0] = _have_loc[1] = 0;
	_max_sift = max_sift <=0 ? 4096 : ((max_sift + 31)/ 32 * 32) ; 
	_pixel_per_sift = 32; //must be 32
	_sift_num_stripe = 1; 
	_sift_per_stripe = 1;
	_sift_per_row = _sift_per_stripe * _sift_num_stripe;
#if defined(CG_SIFTGPU_ENABLED) 
	_use_glsl = use_glsl;
#else
	_use_glsl = 1;
#endif
	_initialized = 0;
}

SiftMatchGL::~SiftMatchGL()
{
	if(s_multiply) delete s_multiply;
	if(s_guided_mult) delete s_guided_mult;
	if(s_col_max) delete s_col_max;
	if(s_row_max) delete s_row_max;
}

void SiftMatchGL::SetMaxSift(int max_sift)
{
	
	max_sift = ((max_sift + 31)/32)*32;
	if(max_sift > GlobalUtil::_texMaxDimGL) max_sift = GlobalUtil::_texMaxDimGL;
	if(max_sift > _max_sift)
	{
		_max_sift = max_sift;
		AllocateSiftMatch();
		_have_loc[0] = _have_loc[1] = 0;
		_id_sift[0] = _id_sift[1] = -1;
		_num_sift[0] = _num_sift[1] = 1;
	}else
	{
		_max_sift = max_sift;
	}

}

void SiftMatchGL::AllocateSiftMatch()
{
	//parameters, number of sift is limited by the texture size
	if(_max_sift > GlobalUtil::_texMaxDimGL) _max_sift = GlobalUtil::_texMaxDimGL;	
	///
	int h = _max_sift / _sift_per_row; 
	int n = (GlobalUtil::_texMaxDimGL + h - 1) / GlobalUtil::_texMaxDimGL; 
	if ( n > 1) {_sift_num_stripe *= n; _sift_per_row *= n; }

	//initialize

	_texDes[0].InitTexture(_sift_per_row * _pixel_per_sift, _max_sift / _sift_per_row, 0,GL_RGBA8);
	_texDes[1].InitTexture(_sift_per_row * _pixel_per_sift, _max_sift / _sift_per_row, 0, GL_RGBA8);
	_texLoc[0].InitTexture(_sift_per_row , _max_sift / _sift_per_row, 0);
	_texLoc[1].InitTexture(_sift_per_row , _max_sift / _sift_per_row, 0);

	if(GlobalUtil::_SupportNVFloat || GlobalUtil::_SupportTextureRG)
	{
		//use single-component texture to save memory
#ifndef GL_R32F
#define GL_R32F 0x822E
#endif
		GLuint format = GlobalUtil::_SupportNVFloat ? GL_FLOAT_R_NV : GL_R32F;
		_texDot.InitTexture(_max_sift, _max_sift, 0, format);
		_texMatch[0].InitTexture(16, _max_sift / 16, 0, format);
		_texMatch[1].InitTexture(16, _max_sift / 16, 0, format);
	}else
	{
		_texDot.InitTexture(_max_sift, _max_sift, 0);
		_texMatch[0].InitTexture(16, _max_sift / 16, 0);
		_texMatch[1].InitTexture(16, _max_sift / 16, 0);
	}

}
void SiftMatchGL::InitSiftMatch()
{
	if(_initialized) return;
	GlobalUtil::InitGLParam();
	if(GlobalUtil::_GoodOpenGL == 0) return;
	AllocateSiftMatch();
#if defined(CG_SIFTGPU_ENABLED) 
	if(_use_glsl)
	{
		LoadSiftMatchShadersGLSL();
	}
	else 
	{
		ProgramCG::InitContext();
		LoadSiftMatchShadersCG();
	}
#else
	LoadSiftMatchShadersGLSL();
#endif
	_initialized = 1; 
}


void SiftMatchGL::SetDescriptors(int index, int num, const unsigned char* descriptors, int id)
{	
	if(_initialized == 0) return;
	if (index > 1) index = 1;
	if (index < 0) index = 0;
	_have_loc[index] = 0;

	//the same feature is already set
	if(id !=-1 && id == _id_sift[index]) return ;
	_id_sift[index] = id;

	if(num > _max_sift) num = _max_sift;

	sift_buffer.resize(num * 128 /4);
	memcpy(&sift_buffer[0], descriptors, 128 * num);
	_num_sift[index] = num; 
	int w = _sift_per_row * _pixel_per_sift;
	int h = (num + _sift_per_row  - 1)/ _sift_per_row; 
	sift_buffer.resize(w * h * 4, 0);
	_texDes[index].SetImageSize(w , h);
	_texDes[index].BindTex(); 
	if(_sift_num_stripe == 1)
	{
		glTexSubImage2D(GlobalUtil::_texTarget, 0, 0, 0, w, h, GL_RGBA,  GL_UNSIGNED_BYTE, &sift_buffer[0]);
	}else
	{
		for(int i = 0; i < _sift_num_stripe; ++i)
		{
			int ws = _sift_per_stripe * _pixel_per_sift;
			int x = i * ws;
			int pos = i * ws * h * 4; 
			glTexSubImage2D(GlobalUtil::_texTarget, 0, x, 0, ws, h, GL_RGBA, GL_UNSIGNED_BYTE, &sift_buffer[pos]);
		}
	}
	_texDes[index].UnbindTex();

}

void SiftMatchGL::SetFeautreLocation(int index, const float* locations, int gap)
{
	if(_num_sift[index] <=0) return;
	int w = _sift_per_row ;
	int h = (_num_sift[index] + _sift_per_row  - 1)/ _sift_per_row; 
	sift_buffer.resize(_num_sift[index] * 2);
	if(gap == 0)
	{
		memcpy(&sift_buffer[0], locations, _num_sift[index] * 2 * sizeof(float));
	}else
	{
		for(int i = 0; i < _num_sift[index]; ++i)
		{
			sift_buffer[i*2] = *locations++;
			sift_buffer[i*2+1]= *locations ++;
			locations += gap;
		}
	}
	sift_buffer.resize(w * h * 2, 0);
	_texLoc[index].SetImageSize(w , h);
	_texLoc[index].BindTex(); 
	if(_sift_num_stripe == 1)
	{
		glTexSubImage2D(GlobalUtil::_texTarget, 0, 0, 0, w, h, GL_LUMINANCE_ALPHA , GL_FLOAT , &sift_buffer[0]);
	}else
	{
		for(int i = 0; i < _sift_num_stripe; ++i)
		{
			int ws = _sift_per_stripe;
			int x = i * ws;
			int pos = i * ws * h * 2; 
			glTexSubImage2D(GlobalUtil::_texTarget, 0, x, 0, ws, h, GL_LUMINANCE_ALPHA , GL_FLOAT, &sift_buffer[pos]);
		}
	}
	_texLoc[index].UnbindTex();
	_have_loc[index] = 1;
}

void SiftMatchGL::SetDescriptors(int index, int num, const float* descriptors, int id)
{	
	if(_initialized == 0) return;
	if (index > 1) index = 1;
	if (index < 0) index = 0;
	_have_loc[index] = 0;

	//the same feature is already set
	if(id !=-1 && id == _id_sift[index]) return ;
	_id_sift[index] = id; 

	if(num > _max_sift) num = _max_sift;

	sift_buffer.resize(num * 128 /4);
	unsigned char * pub = (unsigned char*) &sift_buffer[0];
	for(int i = 0; i < 128 * num; ++i)
	{
		pub[i] = int(512 * descriptors[i] + 0.5);
	}
	_num_sift[index] = num; 
	int w = _sift_per_row * _pixel_per_sift;
	int h = (num + _sift_per_row  - 1)/ _sift_per_row; 
	sift_buffer.resize(w * h * 4, 0);
	_texDes[index].SetImageSize(w, h);
	_texDes[index].BindTex();
	if(_sift_num_stripe == 1)
	{
		glTexSubImage2D(GlobalUtil::_texTarget, 0, 0, 0, w, h, GL_RGBA,  GL_UNSIGNED_BYTE, &sift_buffer[0]);
	}else
	{
		for(int i = 0; i < _sift_num_stripe; ++i)
		{
			int ws = _sift_per_stripe * _pixel_per_sift;
			int x = i * ws;
			int pos = i * ws * h * 4; 
			glTexSubImage2D(GlobalUtil::_texTarget, 0, x, 0, ws, h, GL_RGBA, GL_UNSIGNED_BYTE, &sift_buffer[pos]);
		}
	}
	_texDes[index].UnbindTex();
}

void SiftMatchGL::LoadSiftMatchShadersCG()
{
#if defined(CG_SIFTGPU_ENABLED) 
	ProgramCG * program;
	char buffer[10240];
	ostrstream out(buffer, 10240);
	out <<	"#define SIFT_PER_STRIPE " << _sift_per_stripe << "\n" 
			"#define PIXEL_PER_SIFT " << _pixel_per_sift << "\n"
			"void main(uniform samplerRECT tex1,		\n"
			"uniform samplerRECT tex2,			\n"
			"uniform float2	size,  \n"
			"in		float2	pos : WPOS,	\n"
			"out	float4  result:COLOR0)		\n"
		    "{\n"
		<<	"   float4 val = float4(0, 0, 0, 0), data1, buf;\n"
			"   float2 index = pos.yx; \n"
			"   float2 stripe_size = size.xy * SIFT_PER_STRIPE;\n"
			"   float2 temp_div1 = index / stripe_size;\n"
			"   float2 stripe_index = floor(temp_div1);\n"
			"   index = floor(stripe_size * (temp_div1 - stripe_index));\n"
			"   float2 temp_div2 = index / SIFT_PER_STRIPE;\n"
			"   float2 temp_floor2 = floor(temp_div2);\n"
			"   float2 index_v = temp_floor2 + 0.5;\n "
			"   float2 index_h = SIFT_PER_STRIPE* (temp_div2 - temp_floor2);\n"
			"   float2 tx = (index_h + stripe_index * SIFT_PER_STRIPE)* PIXEL_PER_SIFT + 0.5;\n"
			"   float2 tpos1, tpos2; \n"
			"   float4 tpos = float4(tx, index_v);\n"
			//////////////////////////////////////////////////////
			"   for(int i = 0; i < PIXEL_PER_SIFT; ++i){\n"
			"       buf = texRECT(tex2, tpos.yw);\n"
			"       data1 = texRECT(tex1, tpos.xz);\n"
			"       val += data1 * buf;\n"
			"       tpos.xy = tpos.xy + float2(1.0, 1.0);\n"
			"   }\n"
			"   const float factor = 0.248050689697265625; \n"
			"   result =float4(dot(val, factor.xxxx), index,  0);\n"
			"}"
		<<	'\0';

	s_multiply = program= new ProgramCG(buffer); 

	_param_multiply_tex1 = cgGetNamedParameter(*program, "tex1");
	_param_multiply_tex2 = cgGetNamedParameter(*program, "tex2");
	_param_multiply_size = cgGetNamedParameter(*program, "size");

	out.seekp(ios::beg);
	out <<	"#define SIFT_PER_STRIPE " << _sift_per_stripe << "\n" 
			"#define PIXEL_PER_SIFT " << _pixel_per_sift << "\n"
			"void main(uniform samplerRECT tex1,		\n"
			"uniform samplerRECT tex2,			\n"
			"uniform samplerRECT texL1, \n"
			"uniform samplerRECT texL2, \n"
			"uniform float3x3 H, \n"
			"uniform float3x3 F, \n"
			"uniform float4	size,  \n"
			"in		float2	pos : WPOS,	\n"
			"out	float4  result:COLOR0)		\n"
		    "{\n"
		<<	"   float4 val = float4(0, 0, 0, 0), data1, buf;\n"
			"   float2 index = pos.yx; \n"
			//"   float2 stripe_size = size.xy * SIFT_PER_STRIPE;\n"
			//"   float2 stripe_index = floor(index /stripe_size);\n"
			//"   index = floor(fmod(index,  stripe_size));\n"
			//"   float2 index_v = floor(index / SIFT_PER_STRIPE) + 0.5;\n "
			//"   float2 index_h = fmod(index, SIFT_PER_STRIPE);\n"
			"   float2 stripe_size = size.xy * SIFT_PER_STRIPE;\n"
			"   float2 temp_div1 = index / stripe_size;\n"
			"   float2 stripe_index = floor(temp_div1);\n"
			"   index = floor(stripe_size * (temp_div1 - stripe_index));\n"
			"   float2 temp_div2 = index / SIFT_PER_STRIPE;\n"
			"   float2 temp_floor2 = floor(temp_div2);\n"
			"   float2 index_v = temp_floor2 + 0.5;\n "
			"   float2 index_h = SIFT_PER_STRIPE* (temp_div2 - temp_floor2);\n"
			//read feature location data
			"   float4 tlpos = float4((index_h + stripe_index * SIFT_PER_STRIPE) + 0.5, index_v);\n"
			"   float3 loc1 = float3(texRECT(texL1, tlpos.xz).xw, 1);\n"
			"   float3 loc2 = float3(texRECT(texL2, tlpos.yw).xw, 1);\n"
			//check the guiding homography
			"   float3 hxloc1 = mul(H, loc1);\n"
			"   float2 diff = abs(loc2.xy- hxloc1.xy/hxloc1.z);\n"
			"   float disth = max(diff.x, diff.y);\n"
			"   if(disth > size.z ) {result = float4(0, index, 0); return;}\n"

			//check the guiding fundamental 
			"   float3 fx1 = mul(F, loc1), ftx2 = mul(loc2, F);\n"
			"   float x2tfx1 = dot(loc2, fx1);\n"
			"   float4 temp = float4(fx1.xy, ftx2.xy); \n"
			"   float sampson_error = (x2tfx1 * x2tfx1) / dot(temp, temp);\n"
			"   if(sampson_error > size.w) {result = float4(0, index, 0); return;}\n"

			//compare feature descriptor
			"   float2 tx = (index_h + stripe_index * SIFT_PER_STRIPE)* PIXEL_PER_SIFT + 0.5;\n"
			"   float2 tpos1, tpos2; \n"
			"   float4 tpos = float4(tx, index_v);\n"
			"   for(int i = 0; i < PIXEL_PER_SIFT; ++i){\n"
			"       buf = texRECT(tex2, tpos.yw);\n"
			"       data1 = texRECT(tex1, tpos.xz);\n"
			"       val += data1 * buf;\n"
			"       tpos.xy = tpos.xy + float2(1.0, 1.0);\n"
			"   }\n"
			"   const float factor = 0.248050689697265625; \n"
			"   result =float4(dot(val, factor.xxxx), index,  0);\n"
			"}"
		<<	'\0';

	s_guided_mult = program= new ProgramCG(buffer);

	_param_guided_mult_tex1 = cgGetNamedParameter(*program, "tex1");
	_param_guided_mult_tex2= cgGetNamedParameter(*program, "tex2");
	_param_guided_mult_texl1 = cgGetNamedParameter(*program, "texL1");
	_param_guided_mult_texl2 = cgGetNamedParameter(*program, "texL2");
	_param_guided_mult_h = cgGetNamedParameter(*program, "H");
	_param_guided_mult_f = cgGetNamedParameter(*program, "F");
	_param_guided_mult_param = cgGetNamedParameter(*program, "size");

	//row max
	out.seekp(ios::beg);
	out <<	"#define BLOCK_WIDTH "  << 16 << "\n"
			"void main (uniform samplerRECT tex, \n"
			"uniform float3 param, in float2 pos : WPOS, \n"
			"out float4 result: COLOR0)\n"
			"{\n"
			"   float index = pos.x + floor(pos.y) * BLOCK_WIDTH; \n"
			"   float2 bestv = -1; float imax = -1;\n"
			"   for(float i = 0; i < param.x; i ++){\n "
			"		float v = texRECT(tex, float2(i + 0.5, index)).r; \n"
			"		imax = v > bestv.r ? i : imax; \n "
			"		bestv  = v > bestv.r? float2(v, bestv.r) : max(bestv, v.xx);\n "
			"	}\n"
			"	bestv = acos(min(bestv, 1.0));\n"
			"	if(bestv.x >= param.y || bestv.x >= param.z * bestv.y) imax = -1;\n"
			"	result = float4(imax, bestv, index);\n"
			"}"
		<<  '\0';
	s_row_max = program= new ProgramCG(buffer); 
	_param_rowmax_param = cgGetNamedParameter(*program, "param");

	out.seekp(ios::beg);
	out <<	"#define BLOCK_WIDTH "  << 16 << "\n"
			"void main (uniform samplerRECT tex, \n"
			"uniform float3 param, in float2 pos : WPOS, \n"
			"out float4 result: COLOR0)\n"
			"{\n"
			"	float index = pos.x + floor(pos.y) * BLOCK_WIDTH; \n"
			"	float2 bestv = -1; float imax = -1;\n"
			"	for(float i = 0; i < param.x; i ++){\n "
			"		float v = texRECT(tex, float2(index, i + 0.5)).r; \n"
			"		imax = (v > bestv.r)? i : imax; \n "
			"		bestv  = v > bestv.r? float2(v, bestv.r) : max(bestv, v.xx);\n "
			"	}\n"
			"	bestv = acos(min(bestv, 1.0));\n"
			"	if(bestv.x >= param.y || bestv.x >= param.z * bestv.y) imax = -1;\n"
			"	result = float4(imax, bestv, index);\n"
			"}"
		<<  '\0';
	s_col_max = program =new ProgramCG(buffer); 
	_param_colmax_param = cgGetNamedParameter(*program, "param");
#endif
}

void SiftMatchGL::LoadSiftMatchShadersGLSL()
{
	ProgramGLSL * program;
	char buffer[10240];
	ostrstream out(buffer, 10240);
	out <<  "#pragma optionNV(ifcvt none)\n"
			"#pragma optionNV(unroll all)\n"
	        "#define SIFT_PER_STRIPE " << _sift_per_stripe << ".0\n" 
			"#define PIXEL_PER_SIFT " << _pixel_per_sift << "\n"
			"uniform sampler2DRect tex1, tex2; uniform vec2	size;\n"
			"void main()		\n"
		    "{\n"
		<<	"   vec4 val = vec4(0, 0, 0, 0), data1, buf;\n"
			"   vec2 index = gl_FragCoord.yx; \n"
			"   vec2 stripe_size = size.xy * SIFT_PER_STRIPE;\n"
			"	vec2 temp_div1 = index / stripe_size;\n"
			"   vec2 stripe_index = floor(temp_div1);\n"
			"   index = floor(stripe_size * (temp_div1 - stripe_index));\n"
			"	vec2 temp_div2 = index * vec2(1.0 / float(SIFT_PER_STRIPE));\n"
			"	vec2 temp_floor2 = floor(temp_div2);\n"
			"   vec2 index_v = temp_floor2 + vec2(0.5);\n "
			"   vec2 index_h = vec2(SIFT_PER_STRIPE)* (temp_div2 - temp_floor2);\n"
			"   vec2 tx = (index_h + stripe_index * vec2(SIFT_PER_STRIPE))* vec2(PIXEL_PER_SIFT) + 0.5;\n"
			"   vec2 tpos1, tpos2; \n"
			"	vec4 tpos = vec4(tx, index_v);\n"
			//////////////////////////////////////////////////////
			"   for(int i = 0; i < PIXEL_PER_SIFT; ++i){\n"
			"		buf = texture2DRect(tex2, tpos.yw);\n"
			"		data1 = texture2DRect(tex1, tpos.xz);\n"
			"		val += (data1 * buf);\n"
			"		tpos.xy = tpos.xy + vec2(1.0, 1.0);\n"
			"	}\n"
			"	const float factor = 0.248050689697265625; \n"
			"	gl_FragColor =vec4(dot(val, vec4(factor)), index,  0);\n"
			"}"
		<<	'\0';

	s_multiply = program= new ProgramGLSL(buffer); 

	_param_multiply_tex1 = glGetUniformLocation(*program, "tex1");
	_param_multiply_tex2 = glGetUniformLocation(*program, "tex2");
	_param_multiply_size = glGetUniformLocation(*program, "size");

	out.seekp(ios::beg);
	out <<  "#pragma optionNV(ifcvt none)\n"
			"#pragma optionNV(unroll all)\n"
	        "#define SIFT_PER_STRIPE " << _sift_per_stripe << ".0\n" 
			"#define PIXEL_PER_SIFT " << _pixel_per_sift << "\n"
			"uniform sampler2DRect tex1, tex2;\n"
			"uniform sampler2DRect texL1;\n"
			"uniform sampler2DRect texL2; \n"
			"uniform mat3 H; \n"
			"uniform mat3 F; \n"
			"uniform vec4	size; \n"
			"void main()		\n"
		    "{\n"
		<<	"   vec4 val = vec4(0, 0, 0, 0), data1, buf;\n"
			"   vec2 index = gl_FragCoord.yx; \n"
			"   vec2 stripe_size = size.xy * SIFT_PER_STRIPE;\n"
			"	vec2 temp_div1 = index / stripe_size;\n"
			"   vec2 stripe_index = floor(temp_div1);\n"
			"   index = floor(stripe_size * (temp_div1 - stripe_index));\n"
			"	vec2 temp_div2 = index  * vec2(1.0/ float(SIFT_PER_STRIPE));\n"
			"	vec2 temp_floor2 = floor(temp_div2);\n"
			"   vec2 index_v = temp_floor2 + vec2(0.5);\n "
			"   vec2 index_h = vec2(SIFT_PER_STRIPE)* (temp_div2 - temp_floor2);\n"
			
			//read feature location data
			"   vec4 tlpos = vec4((index_h + stripe_index * vec2(SIFT_PER_STRIPE)) + 0.5, index_v);\n"
			"   vec3 loc1 = vec3(texture2DRect(texL1, tlpos.xz).xw, 1.0);\n"
			"   vec3 loc2 = vec3(texture2DRect(texL2, tlpos.yw).xw, 1.0);\n"
			
			//check the guiding homography
			"   vec3 hxloc1 = H* loc1;\n"
			"   vec2 diff = abs(loc2.xy- (hxloc1.xy/hxloc1.z));\n"
			"   float disth = max(diff.x, diff.y);\n"
			"   if(disth > size.z ) {gl_FragColor = vec4(0, index, 0); return;}\n"

			//check the guiding fundamental 
			"   vec3 fx1 = (F * loc1), ftx2 = (loc2 * F);\n"
			"   float x2tfx1 = dot(loc2, fx1);\n"
			"   vec4 temp = vec4(fx1.xy, ftx2.xy); \n"
			"   float sampson_error = (x2tfx1 * x2tfx1) / dot(temp, temp);\n"
			"   if(sampson_error > size.w) {gl_FragColor = vec4(0, index, 0); return;}\n"

			//compare feature descriptor
			"   vec2 tx = (index_h + stripe_index * SIFT_PER_STRIPE)* vec2(PIXEL_PER_SIFT) + 0.5;\n"
			"   vec2 tpos1, tpos2; \n"
			"	vec4 tpos = vec4(tx, index_v);\n"
			"   for(int i = 0; i < PIXEL_PER_SIFT; ++i){\n"
			"		buf = texture2DRect(tex2, tpos.yw);\n"
			"		data1 = texture2DRect(tex1, tpos.xz);\n"
			"		val += data1 * buf;\n"
			"		tpos.xy = tpos.xy + vec2(1.0, 1.0);\n"
			"	}\n"
			"	const float factor = 0.248050689697265625; \n"
			"	gl_FragColor =vec4(dot(val, vec4(factor)), index,  0);\n"
			"}"
		<<	'\0';

	s_guided_mult = program= new ProgramGLSL(buffer);

	_param_guided_mult_tex1 = glGetUniformLocation(*program, "tex1");
	_param_guided_mult_tex2= glGetUniformLocation(*program, "tex2");
	_param_guided_mult_texl1 = glGetUniformLocation(*program, "texL1");
	_param_guided_mult_texl2 = glGetUniformLocation(*program, "texL2");
	_param_guided_mult_h = glGetUniformLocation(*program, "H");
	_param_guided_mult_f = glGetUniformLocation(*program, "F");
	_param_guided_mult_param = glGetUniformLocation(*program, "size");

	//row max
	out.seekp(ios::beg);
	out <<	"#define BLOCK_WIDTH 16.0\n"
			"uniform sampler2DRect tex;	uniform vec3 param;\n"
			"void main ()\n"
			"{\n"
			"	float index = gl_FragCoord.x + floor(gl_FragCoord.y) * BLOCK_WIDTH; \n"
			"	vec2 bestv = vec2(-1.0); float imax = -1.0;\n"
			"	for(float i = 0.0; i < param.x; i ++){\n "
			"		float v = texture2DRect(tex, vec2(i + 0.5, index)).r; \n"
			"		imax = v > bestv.r ? i : imax; \n "
			"		bestv  = v > bestv.r? vec2(v, bestv.r) : max(bestv, vec2(v));\n "
			"	}\n"
			"	bestv = acos(min(bestv, 1.0));\n"
			"	if(bestv.x >= param.y || bestv.x >= param.z * bestv.y) imax = -1.0;\n"
			"	gl_FragColor = vec4(imax, bestv, index);\n"
			"}"
		<<  '\0';
	s_row_max = program= new ProgramGLSL(buffer); 
	_param_rowmax_param = glGetUniformLocation(*program, "param");

	out.seekp(ios::beg);
	out <<	"#define BLOCK_WIDTH 16.0\n"
			"uniform sampler2DRect tex; uniform vec3 param;\n"
			"void main ()\n"
			"{\n"
			"	float index = gl_FragCoord.x + floor(gl_FragCoord.y) * BLOCK_WIDTH; \n"
			"	vec2 bestv = vec2(-1.0); float imax = -1.0;\n"
			"	for(float i = 0.0; i < param.x; i ++){\n "
			"		float v = texture2DRect(tex, vec2(index, i + 0.5)).r; \n"
			"		imax = (v > bestv.r)? i : imax; \n "
			"		bestv  = v > bestv.r? vec2(v, bestv.r) : max(bestv, vec2(v));\n "
			"	}\n"
			"	bestv = acos(min(bestv, 1.0));\n"
			"	if(bestv.x >= param.y || bestv.x >= param.z * bestv.y) imax = -1.0;\n"
			"	gl_FragColor = vec4(imax, bestv, index);\n"
			"}"
		<<  '\0';
	s_col_max = program =new ProgramGLSL(buffer); 
	_param_colmax_param = glGetUniformLocation(*program, "param");


}

int  SiftMatchGL::GetGuidedSiftMatch(int max_match, int match_buffer[][2], float H[3][3], float F[3][3],
									 float distmax, float ratiomax, float hdistmax, float fdistmax, int mbm)
{

	int dw = _num_sift[1];
	int dh = _num_sift[0]; 
	if(_initialized ==0) return 0;
	if(dw <= 0 || dh <=0) return 0;
	if(_have_loc[0] == 0 || _have_loc[1] == 0) return 0;

	FrameBufferObject fbo;
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	_texDot.SetImageSize(dw, dh);


	//data
	_texDot.AttachToFBO(0);
	_texDot.FitTexViewPort();
	glActiveTexture(GL_TEXTURE0);
	_texDes[0].BindTex();
	glActiveTexture(GL_TEXTURE1);
	_texDes[1].BindTex();
	glActiveTexture(GL_TEXTURE2);
	_texLoc[0].BindTex();
	glActiveTexture(GL_TEXTURE3);
	_texLoc[1].BindTex();

	//multiply the descriptor matrices
	s_guided_mult->UseProgram();


	//set parameters glsl/cg
	float dot_param[4] = {(float)_texDes[0].GetDrawHeight(), (float) _texDes[1].GetDrawHeight(), hdistmax, fdistmax};
#if defined(CG_SIFTGPU_ENABLED) 
	if(_use_glsl)
	{
#endif
		glUniform1i(_param_guided_mult_tex1, 0);
		glUniform1i(_param_guided_mult_tex2, 1);
		glUniform1i(_param_guided_mult_texl1, 2);
		glUniform1i(_param_guided_mult_texl2, 3);
		glUniformMatrix3fv(_param_guided_mult_h, 1, GL_TRUE, H[0]);
		glUniformMatrix3fv(_param_guided_mult_f, 1, GL_TRUE, F[0]);
		glUniform4fv(_param_guided_mult_param, 1, dot_param);
#if defined(CG_SIFTGPU_ENABLED) 
	}else
	{

		cgGLSetTextureParameter(_param_guided_mult_tex1, _texDes[0]);	
		cgGLEnableTextureParameter(_param_guided_mult_tex1);
		cgGLSetTextureParameter(_param_guided_mult_tex2, _texDes[1]);	
		cgGLEnableTextureParameter(_param_guided_mult_tex2);
		cgGLSetTextureParameter(_param_guided_mult_texl1, _texLoc[0]);	
		cgGLEnableTextureParameter(_param_guided_mult_texl1);
		cgGLSetTextureParameter(_param_guided_mult_texl2, _texLoc[1]);	
		cgGLEnableTextureParameter(_param_guided_mult_texl2);
		cgSetMatrixParameterfr(_param_guided_mult_h, H[0]);
		cgSetMatrixParameterfr(_param_guided_mult_f, F[0]);
		cgGLSetParameter4fv(_param_guided_mult_param, dot_param);
	}
#endif


	_texDot.DrawQuad();

	GLTexImage::UnbindMultiTex(4);

	return GetBestMatch(max_match, match_buffer, distmax, ratiomax, mbm);
}

int SiftMatchGL::GetBestMatch(int max_match, int match_buffer[][2], float distmax, float ratiomax, int mbm)
{

	glActiveTexture(GL_TEXTURE0);
	_texDot.BindTex();

	//readback buffer
	sift_buffer.resize(_num_sift[0] + _num_sift[1] + 16);
	float * buffer1 = &sift_buffer[0], * buffer2 = &sift_buffer[_num_sift[0]];

	//row max
	_texMatch[0].AttachToFBO(0);
	_texMatch[0].SetImageSize(16, ( _num_sift[0] + 15) / 16);
	_texMatch[0].FitTexViewPort();

	///set parameter glsl/cg
	s_row_max->UseProgram();
#if defined(CG_SIFTGPU_ENABLED) 
	if(_use_glsl) glUniform3f(_param_rowmax_param, (float)_num_sift[1], distmax, ratiomax);
	else cgGLSetParameter3f(_param_rowmax_param, (float)_num_sift[1], distmax, ratiomax);
#else
	glUniform3f(_param_rowmax_param, (float)_num_sift[1], distmax, ratiomax);
#endif

	_texMatch[0].DrawQuad();
	glReadPixels(0, 0, 16, (_num_sift[0] + 15)/16, GL_RED, GL_FLOAT, buffer1);

	//col max
	if(mbm)
	{
		_texMatch[1].AttachToFBO(0);
		_texMatch[1].SetImageSize(16, (_num_sift[1] + 15) / 16);
		_texMatch[1].FitTexViewPort();
		//set parameter glsl/cg
		s_col_max->UseProgram();
#if defined(CG_SIFTGPU_ENABLED) 
		if(_use_glsl) glUniform3f(_param_rowmax_param, (float)_num_sift[0], distmax, ratiomax);
		else cgGLSetParameter3f(_param_colmax_param, (float)_num_sift[0], distmax, ratiomax);
#else
		glUniform3f(_param_rowmax_param, (float)_num_sift[0], distmax, ratiomax);
#endif
		_texMatch[1].DrawQuad();
		glReadPixels(0, 0, 16, (_num_sift[1] + 15) / 16, GL_RED, GL_FLOAT, buffer2);
	}


	//unload
#if defined(CG_SIFTGPU_ENABLED) 
	if(_use_glsl)
	{
#endif
		glUseProgram(0);
#if defined(CG_SIFTGPU_ENABLED) 
	}else
	{
		cgGLUnbindProgram(ProgramCG::_FProfile);
		cgGLDisableProfile(ProgramCG::_FProfile);
	}
#endif
	GLTexImage::UnbindMultiTex(2);
	GlobalUtil::CleanupOpenGL();

	//write back the matches
	int nmatch = 0, j ;
	for(int i = 0; i < _num_sift[0] && nmatch < max_match; ++i)
	{
		j = int(buffer1[i]);
		if( j>= 0 && (!mbm ||int(buffer2[j]) == i))
		{
			match_buffer[nmatch][0] = i;
			match_buffer[nmatch][1] = j;
			nmatch++;
		}
	}
	return nmatch;
}

int  SiftMatchGL::GetSiftMatch(int max_match, int match_buffer[][2], float distmax, float ratiomax, int mbm)
{
	int dw = _num_sift[1];
	int dh =  _num_sift[0]; 
	if(_initialized ==0) return 0;
	if(dw <= 0 || dh <=0) return 0;

	FrameBufferObject fbo;
	glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
	_texDot.SetImageSize(dw, dh);

	//data
	_texDot.AttachToFBO(0);
	_texDot.FitTexViewPort();
	glActiveTexture(GL_TEXTURE0);
	_texDes[0].BindTex();
	glActiveTexture(GL_TEXTURE1);
	_texDes[1].BindTex();

	//////////////////	
	//multiply the descriptor matrices
	s_multiply->UseProgram();
	//set parameters
	float heights[2] = {(float)_texDes[0].GetDrawHeight(), (float)_texDes[1].GetDrawHeight()};

#if defined(CG_SIFTGPU_ENABLED) 
	if(_use_glsl)
	{
#endif
		glUniform1i(_param_multiply_tex1, 0);
		glUniform1i(_param_multiply_tex2 , 1);
		glUniform2fv(_param_multiply_size, 1, heights);

#if defined(CG_SIFTGPU_ENABLED) 
	}else
	{
		cgGLSetTextureParameter(_param_multiply_tex1, _texDes[0]);	
		cgGLEnableTextureParameter(_param_multiply_tex1);
		cgGLSetTextureParameter(_param_multiply_tex2, _texDes[1]);	
		cgGLEnableTextureParameter(_param_multiply_tex2);
		cgGLSetParameter2fv(_param_multiply_size, heights);
	}
#endif

	_texDot.DrawQuad();

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GlobalUtil::_texTarget, 0);

	return GetBestMatch(max_match, match_buffer, distmax, ratiomax, mbm);
}


int SiftMatchGPU::_CreateContextGL()
{
	//use GLUT to create an OpenGL Context?
    if (__language >= SIFTMATCH_CUDA && GlobalUtil::_NoOpenGLinCUDA) {}
	else if(!GlobalUtil::CreateWindowEZ()) return 0;
	return VerifyContextGL();
}


int SiftMatchGPU::_VerifyContextGL()
{
	if(__matcher) return GlobalUtil::_GoodOpenGL;
	
#ifdef CUDA_SIFTGPU_ENABLED
    if(((__language == SIFTMATCH_SAME_AS_SIFTGPU && GlobalUtil::_UseCUDA) || __language >= SIFTMATCH_CUDA) 
        && SiftMatchCU::CheckCudaDevice (GlobalUtil::_CudaDeviceIndex))
    {
		__language = SIFTMATCH_CUDA;
		__matcher = new SiftMatchCU(__max_sift);
	}else
#else
    if((__language == SIFTMATCH_SAME_AS_SIFTGPU && GlobalUtil::_UseCUDA) || __language >= SIFTMATCH_CUDA) 
    {
	    std::cerr	<< "---------------------------------------------------------------------------\n"
				    << "CUDA not supported in this binary! To enable it, please use SiftGPU_CUDA_Enable\n" 
				    << "Project for VS2005+ or set siftgpu_enable_cuda to 1 in makefile\n"
				    << "----------------------------------------------------------------------------\n";
    }
#endif
#ifdef CG_SIFTGPU_ENABLED
	if((__language == SIFTMATCH_SAME_AS_SIFTGPU && GlobalUtil::_UseGLSL == 0) || __language == SIFTMATCH_CG)
	{
		__language = SIFTMATCH_CG;
		__matcher = new SiftMatchGL(__max_sift, 0);
	}else
#endif
	{
		__language = SIFTMATCH_GLSL;
		__matcher = new SiftMatchGL(__max_sift, 1);
	}

	if(GlobalUtil::_verbose)	std::cout<<"[SiftMatchGPU Language]: "<<
		(__language == SIFTMATCH_CUDA? "CUDA" : (__language == SIFTMATCH_GLSL?"GLSL" : "CG")) <<"\n\n";

	__matcher->InitSiftMatch();
	return GlobalUtil::_GoodOpenGL;
}

void* SiftMatchGPU::operator new (size_t  size){
  void * p = malloc(size);
  if (p == 0)  
  {
	  const std::bad_alloc ba;
	  throw ba; 
  }
  return p; 
}


SiftMatchGPU::SiftMatchGPU(int max_sift)
{
	__max_sift = max(max_sift, 1024);
	__language = 0;
	__matcher = NULL;
}

void SiftMatchGPU::SetLanguage(int gpu_language)
{
	if(__matcher) return;
    ////////////////////////

	__language = gpu_language;

#ifdef CUDA_SIFTGPU_ENABLED
	if(__language >= SIFTMATCH_CUDA)
    {
        GlobalUtil::_CudaDeviceIndex = __language - SIFTMATCH_CUDA; 
    }
#else
    if(__language >= SIFTMATCH_CUDA) __language = SIFTMATCH_GLSL; 
#endif
}

void SiftMatchGPU::SetDeviceParam(int argc, char**argv)
{
    if(__matcher) return;

    GlobalUtil::SetDeviceParam(argc, argv);
}

void SiftMatchGPU::SetMaxSift(int max_sift)
{
	if(__matcher)	__matcher->SetMaxSift(max(128, max_sift));
	else __max_sift = max(128, max_sift);
}

SiftMatchGPU::~SiftMatchGPU()
{
	if(__matcher) delete __matcher;
}

void SiftMatchGPU::SetDescriptors(int index, int num, const unsigned char* descriptors, int id)
{
	__matcher->SetDescriptors(index, num,  descriptors, id);
}

void SiftMatchGPU::SetDescriptors(int index, int num, const float* descriptors, int id)
{
	__matcher->SetDescriptors(index, num, descriptors, id);
}

void SiftMatchGPU::SetFeautreLocation(int index, const float* locations, int gap)
{
	__matcher->SetFeautreLocation(index, locations, gap);

}
int  SiftMatchGPU::GetGuidedSiftMatch(int max_match, int match_buffer[][2], float H[3][3], float F[3][3], 
				float distmax, float ratiomax, float hdistmax, float fdistmax, int mutual_best_match)
{
	if(H == NULL && F == NULL)
	{
		return __matcher->GetSiftMatch(max_match, match_buffer, distmax, ratiomax, mutual_best_match);
	}else
	{
		float Z[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, ti = (1.0e+20F);

		return __matcher->GetGuidedSiftMatch(max_match, match_buffer, H? H : Z, F? F : Z,
			distmax, ratiomax, H? hdistmax: ti,  F? fdistmax: ti, mutual_best_match);
	}
}

int  SiftMatchGPU::GetSiftMatch(int max_match, int match_buffer[][2], float distmax, float ratiomax, int mutual_best_match)
{
	return __matcher->GetSiftMatch(max_match, match_buffer, distmax, ratiomax, mutual_best_match);
}

SiftMatchGPU* CreateNewSiftMatchGPU(int max_sift)
{
	return new SiftMatchGPU(max_sift);
}

