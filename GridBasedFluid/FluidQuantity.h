#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

// For definition of spatial coordinate and array index, you can find 
// my note for reference.
class FluidQuantity
{
public:
	// The sequence of value.
	vector<double> mSrcVec;
	vector<double> mDstVec;

	// Size of this DataStructure in 2D.
	// I define that every pixel is a cell in grid.
	// e.g. For a 400 x 400 image, there are 400 x 400 density quantities.
	// Meanwhile, there are 401 x 400 u quantities and 400 x 401 v quantities.
	int mWidth;
	int mHeight;

	// Index offset:
	// e.g. centered quantities: (0.5, 0.5); jittered quantities: (0.0, 0.5) or (0.5, 0.0)
	double mOffsetX;
	double mOffsetY;

	// A grid cell size:
	// Default: it should be size of a pixel.
	double mGridCellSize;

	FluidQuantity(){}

	FluidQuantity(int w, int h, double ox, double oy, double hx)
		: mWidth(w), mHeight(h), mOffsetX(ox), mOffsetY(oy), mGridCellSize(hx)
	{
		mSrcVec.resize(mWidth * mHeight);
		mDstVec.resize(mWidth * mHeight);
	}

	double at(int x, int y) const {
		return mSrcVec.at(x + y * mWidth);
	}

	double &at(int x, int y) {
		return mSrcVec.at(x + y * mWidth);
	}
	// The input x,y is a spatial coordinate.
	// The output is the interploated value at [x, y].
	double BilinearInterpolation(double x, double y) const {
		// Benediket Bitterli's code:
		// It is correct by my test.
		// However, I have no idea why it is correct.
		x = std::min(std::max(x - mOffsetX, 0.0), mWidth - 1.001);
		y = std::min(std::max(y - mOffsetY, 0.0), mHeight - 1.001);
		int ix = (int)x;
		int iy = (int)y;
		x -= ix;
		y -= iy;

		double x00 = this->at(ix + 0, iy + 0);
		double x10 = this->at(ix + 1, iy + 0);
		double x01 = this->at(ix + 0, iy + 1);
		double x11 = this->at(ix + 1, iy + 1);

		return lerp(lerp(x00, x10, x), lerp(x01, x11, x), y);
	}

	// Advect grid in velocity field u, v with given timestep.
	void advect(double timestep, const FluidQuantity &u, const FluidQuantity &v)
	{
		for (int iy = 0, idx = 0; iy < mHeight; iy++)
		{
			for (int ix = 0; ix < mWidth; ix++)
			{
				// Spatial position of data at [x, y].
				double x = ix + mOffsetX;
				double y = iy + mOffsetY;

				// Semi-Lagrangian view:
				// Get particle's position [x', y'].
				euler(x, y, timestep, u, v);

				// Get the quantity value at [x', y'].
				// Use it as next frame's quantity value at [x, y].
				mDstVec.at(idx) = this->BilinearInterpolation(x, y);
			}
		}
	}
	
	void flip() {
		swap(mSrcVec, mDstVec);
	}

	// Sets fluid quantity inside the given rect to the value
	void addInflow(double x0, double y0, double x1, double y1, double val)
	{
		int ix0 = (int)(x0 - mOffsetX);
		int iy0 = (int)(y0 - mOffsetY);
		int ix1 = (int)(x1 - mOffsetX);
		int iy1 = (int)(y1 - mOffsetY);

		for (int y = max(iy0, 0); y < min(iy1, mHeight); y++)
		{
			for (int x = max(ix0, 0); x < min(ix1, mWidth); x++)
			{
				if (fabs(mSrcVec.at(x + y * mWidth)) < fabs(val))
				{
					mSrcVec.at(x + y * mWidth) = val;
				}
			}
		}
	}
	
	// Semi-Lagrangian Euler method:
	// According to the velocity field 'u', 'v' and spatial position [x, y], we are going to find
	// [x', y'] at last frame.
	void euler(double& x, double& y, double timestep, const FluidQuantity& u, const FluidQuantity& v)
	{
		// Find out the velocity at [x, y].
		double uVel = u.BilinearInterpolation(x, y);
		double vVel = v.BilinearInterpolation(x, y);

		// Calculate the [x', y'].
		x -= uVel * timestep;
		y -= vVel * timestep;
	}

	// This class represents a type of data in the MAC grid or on the grid.
	// e.g. all 'u' values or all 'density' values
	double lerp(double a, double b, double x) const {
		return a * (1.0 - x) + b * x;
	}
};


