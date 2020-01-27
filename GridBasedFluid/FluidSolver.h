#pragma once
#include "FluidQuantity.h"
#include <iostream>
// Sets up the fluid quantities;
// Forces incompressibility;
// Performs advection;
// Adds inflows;
class FluidSolver
{
public:
	FluidQuantity mDensity; // We can imagine that this is the density of ink.
	FluidQuantity mU;
	FluidQuantity mV;
	vector<double> mPressure;
	vector<double> mRhs;

	int mWidth;
	int mHeight;

	double mGridCellSize;
	double mFluidDensity; // This is the density of the fluid.

	FluidSolver(int w, int h, double density)
		: mWidth(w), mHeight(h), mFluidDensity(density)
	{
		mDensity = FluidQuantity(mWidth, mHeight, 0.5, 0.5, 1.0);
		mU = FluidQuantity(mWidth + 1, mHeight, 0.0, 0.5, 1.0);
		mV = FluidQuantity(mWidth, mHeight + 1, 0.5, 0.0, 1.0);

		mPressure.resize(mWidth * mHeight);
		mRhs.resize(mWidth * mHeight);
	}

	void addInflow(double x, double y, double w, double h, double dVal, double uVal, double vVal)
	{
		mDensity.addInflow(x, y, x + w, y + h, dVal);
		mU.addInflow(x, y, x + w, y + h, uVal);
		mV.addInflow(x, y, x + w, y + h, vVal);
	}

	double maxTimestep() {
		double maxVelocity = 0.0;
		for (int y = 0; y < mHeight; y++)
		{
			for (int x = 0; x < mWidth; x++)
			{
				double u = mU.BilinearInterpolation(x + 0.5, y + 0.5);
				double v = mV.BilinearInterpolation(x + 0.5, y + 0.5);
				double velocity = sqrt(u * u + v * v);
				maxVelocity = max(maxVelocity, velocity);
			}
		}
		double maxTimestep = 2.0 / maxVelocity;

		return min(maxTimestep, 1.0);
	}

	// This is not the method to create Rhs in the book.
	// Rhs = negative divergence.
	// A * P = Rhs.
	void buildRhs() {
		// double scale = 1.0;
		mRhs.resize(mWidth * mHeight);
		for (int y = 0, idx = 0; y < mHeight; y++)
		{
			for (int x = 0; x < mWidth; x++, idx++)
			{
				mRhs.at(idx) = -(mU.at(x + 1, y) - mU.at(x, y) + mV.at(x, y + 1) - mV.at(x, y));
			}
		}
	}

	void project(int limit, double timestep)
	{
		double scale = timestep / (mFluidDensity);

		double maxDelta;
		for (int iter = 0; iter < limit; iter++) {
			maxDelta = 0.0;
			for (int y = 0, idx = 0; y < mHeight; y++) {
				for (int x = 0; x < mWidth; x++, idx++) {
					int idx = x + y * mWidth;

					double diag = 0.0, offDiag = 0.0;

					/* Here we build the matrix implicitly as the five-point
					 * stencil. Grid borders are assumed to be solid, i.e.
					 * there is no fluid outside the simulation domain.
					 */
					if (x > 0) {
						diag += scale;
						offDiag -= scale * mPressure[idx - 1];
					}
					if (y > 0) {
						diag += scale;
						offDiag -= scale * mPressure[idx - mWidth];
					}
					if (x < mWidth - 1) {
						diag += scale;
						offDiag -= scale * mPressure[idx + 1];
					}
					if (y < mHeight - 1) {
						diag += scale;
						offDiag -= scale * mPressure[idx + mWidth];
					}

					double newP = (mRhs[idx] - offDiag) / diag;

					maxDelta = max(maxDelta, fabs(mPressure[idx] - newP));

					mPressure[idx] = newP;
				}
			}

			if (maxDelta < 1e-5)
			{
				cout << "Exiting solver after" << iter << "iterations, maximum change is " << maxDelta << endl;
				return;
			}
		}
		cout << "Exceeded budget of " << limit << "iterations, maximum change was " << maxDelta << endl;
	}

	void applyPressure(double timestep)
	{
		double scale = timestep / (mFluidDensity * 1.0);

		// For fluid cells:
		for (int y = 0, idx = 0; y < mHeight; y++)
		{
			for (int x = 0; x < mWidth; x++, idx++)
			{
				mU.at(x, y) -= scale * mPressure.at(idx);
				mU.at(x + 1, y) += scale * mPressure.at(idx);
				mV.at(x, y) -= scale * mPressure.at(idx);
				mV.at(x, y + 1) += scale * mPressure.at(idx);
			}
		}

		// For solid boundary, set velocity values as 0.0;
		for (int y = 0; y < mHeight; y++)
		{
			mU.at(0, y) = 0.0;
			mU.at(mWidth, y) = 0.0;
		}

		for (int x = 0; x < mWidth; x++)
		{
			mV.at(x, 0) = 0.0;
			mV.at(x, mHeight) = 0.0;
		}
	}

	void update(double timestep)
	{
		buildRhs();
		project(600, timestep);
		applyPressure(timestep);

		mDensity.advect(timestep, mU, mV);
		mU.advect(timestep, mU, mV);
		mV.advect(timestep, mU, mV);

		mDensity.flip();
		mU.flip();
		mV.flip();
	}

	void toImage(unsigned char *rgb) {
		for (int i = 0; i < mWidth * mHeight; i++)
		{
			// The larger the density, the darker the area.
			int greyScale = (int)((1.0 - mDensity.mSrcVec.at(i)) * 255.0);
			greyScale = max(min(greyScale, 255), 0);
			rgb[i * 3 + 0] = greyScale;
			rgb[i * 3 + 1] = greyScale;
			rgb[i * 3 + 2] = greyScale;
		}
	}
};

