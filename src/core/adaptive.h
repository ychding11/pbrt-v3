
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CORE_ADAPTIVE_H
#define PBRT_CORE_ADAPTIVE_H

// core/adaptive.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include <set>
#include <iterator>
#include <vector>

#define FOR(i,a,b) for(int i=a;i<=b;i++)
#define mp make_pair
#define fst first
#define snd second
#define pb push_back

#define MULT 1

namespace pbrt {

	using namespace std;

struct CDFAdaptive
{
	typedef Vector3<Float> Point;
	typedef Vector3<Float> Vector;

	CDFAdaptive(const float*cdf, int count, int type, int mSize, float mDist, float mRDist);
	void printfCdfAdap(void);
	void toArrays(std::vector<float>& cdfVal, std::vector<int>& cdfIndex);
	void updateFunc(const std::vector<float>&f, std::vector<float>&func);
		
private:
			int trianglArea(int i, int j, int k, const float *f);
			void DouglasPeucker(const float *f, int n);
			void RadialDistance(const float *f, int n);
			void ReumannWitkam(const float *f, int n);
			void Opheim(const float *f, int n);
			void Lang(const float *f, int n);
			void VisvalingamWhyatt(const float *f, int n);
			void PerpendicularDistance(const float *f, int n);
public:
			int count;
			int maxSize;
			float minPDist;
			float minRDist;
			set< pair<int, float> > cdf;

			float distance(int i, int j, const float *f)
			{
				Point a = Point(float(i), f[i], 1.0);
				Point b = Point(float(j), f[j], 1.0);
				Vector dist = a - b;
				return floor(dist.Length() * MULT);
			}
};

	// Monte Carlo Utility Declarations
	struct Distribution1DAdaptive
	{
		// Distribution1D Public Methods
		Distribution1DAdaptive(const float *f, int n, int type, int mSize, float mDist, float mRDist);
		Distribution1DAdaptive(Distribution1D *dist, int tp, int mSize, float mDist, float mRDist);

		~Distribution1DAdaptive() {}
		float SampleContinuous(float u, float *pdf, int *off);
		float Pdf(float u, int &adaptiveU);

		int uniformCount;
		int adaptiveCount;
		float funcInt;
		std::vector<Float> func;
		std::vector<Float> cdf;
		std::vector<int> index;
	};

	struct Distribution2DAdaptive
	{
		// Distribution2DAdaptive Public Methods
		Distribution2DAdaptive(const float *data, int nu, int nv, int tp, int mSize, float mDist, float mRDist);
		virtual ~Distribution2DAdaptive() {};
		void SampleContinuous(float u0, float u1, float uv[2], float *pdf);
		float Pdf(float u, float v);
		int getMargCount() const { return marginalCount; };
		int getCondCount() const { return conditionalCount; };

	private:
		// Distribution2DAdaptive Private Data
		vector<std::unique_ptr<Distribution1DAdaptive> > pConditionalV;
		std::unique_ptr<Distribution1DAdaptive> pMarginal;
		int marginalCount;
		int conditionalCount;
	};

	inline float calcDistance(std::pair<int, float>p1, std::pair<int, float>p2, std::pair<int, float>p)
	{
		float dist = fabs((p2.fst - p1.fst)*(p1.snd - p.snd)
			- (p1.fst - p.fst)*(p2.snd - p1.snd));
		return dist / sqrt((p2.fst - p1.fst)*(p2.fst - p1.fst)
			+ (p2.snd - p1.snd)*(p2.snd - p1.snd));
	}

} //namespace

#endif // PBRT_CORE_ADAPTIVE_H
