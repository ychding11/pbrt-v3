
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

#ifndef PBRT_MATERIALS_MEASURED_ADAPTIVE
#define PBRT_MATERIALS_MEASURED_ADAPTIVE

// materials/measuredAdaptive.h*
#include "pbrt.h"
#include "material.h"
#include "reflection.h"
#include "adaptive.h"
//#include "kdtree.h"

#include <vector>
#include <map>

namespace pbrt {

#ifndef M_PI
#define M_PI (3.14159265358939723846f)
#endif
	class AdaptiveHalfangleBRDF : public BxDF
	{
		typedef Vector3<Float> Point;
		typedef Vector3<Float> Vector;

	public:
		AdaptiveHalfangleBRDF(const std::shared_ptr<float> &d,
			uint32_t nth, uint32_t ntd, uint32_t npd,
			uint32_t mto, uint32_t mpo, uint32_t mth,
			uint32_t mph, vector<Distribution2DAdaptive*> dist);

		~AdaptiveHalfangleBRDF()
		{
			//if (brdf) delete[] brdf;
			for (uint32_t i = 0; i < distribution.size(); ++i)
				if (distribution[i]) delete distribution[i];
		}
		Spectrum f(const Vector &wo, const Vector &wi) const override;
		Spectrum Sample_f(const Vector &wo, Vector *wi, const Point2f &sample, float *pdf, BxDFType *sampledType = nullptr) const override;
		float Pdf(const Vector &wi, const Vector &wo) const override;
		std::string ToString() const override;

		const std::shared_ptr<float> brdf;
		uint32_t nThetaH, nThetaD, nPhiD;
		uint32_t mThetaO, mPhiO, mThetaH, mPhiH;
		vector<Distribution2DAdaptive* > distribution;
	};

	// MeasuredAdaptiveMaterial Declarations
	class MeasuredAdaptiveMaterial : public Material
	{
		typedef Vector3<Float> Vector;
	public:
		// MeasuredAdaptiveMaterial Public Methods
		MeasuredAdaptiveMaterial(const string &filename, const std::shared_ptr<Texture<float> >& bump, int type,
			int mSize, float mPDist, float mRDist);
		~MeasuredAdaptiveMaterial()
		{
			//if (regularHalfangleData) delete[]regularHalfangleData;
			for (uint32_t i = 0; i < distribution.size(); ++i)
				if (distribution[i]) delete distribution[i];

		}

		// OVERRIDE
		void ComputeScatteringFunctions(SurfaceInteraction *si, MemoryArena &arena, TransportMode mode,
			bool allowMultipleLobes) const override;
	private:
		void mapToViewHalfangle(float *tmpData, float *finalData);
		int lookup_brdf_val(float *brdf, double thetaOut, double phiOut, double thetaHalf, double phiHalf);
		void loadAndAnalyzeBRDF(const string &filename, int type, int mSize, float mPDist, float mRDist);

	private:
		// MeasuredAdaptiveMaterial Private Data
		const uint32_t nThetaH = 90, nThetaD = 90, nPhiD = 180,
			mThetaO = 32, mPhiO = 16, mThetaH = 256, mPhiH = 32;
		//float *regularHalfangleData;
		std::shared_ptr<float> regularHalfAngleData;
		vector<Distribution2DAdaptive *> distribution;
		std::shared_ptr<Texture<Float>> bumpMap;

	public:
		// Class member
		static std::map<string, std::shared_ptr<float>> sLoadedRegularHalfAngleAdaptive;
		static std::map<string, vector <Distribution2DAdaptive*> > sLoadedDistributionAdaptive;
	};

	MeasuredAdaptiveMaterial* CreateMeasuredAdaptiveMaterial(const TextureParams &mp);
} //namespace

#endif // PBRT_MATERIALS_MEASURED_ADAPTIVE
