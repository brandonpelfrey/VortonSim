#ifndef UniformGrid_H
#define UniformGrid_H

#include "Vec3.h"
#include "Mat3.h"

template<typename T> T zero();
template<> float zero() { return 0.f; }
template<> Vec3 zero() { return Vec3(0,0,0); }
template<> Mat3 zero() { Vec3 z(0,0,0); return Mat3(z,z,z); }

template<typename T>
class UniformGrid {
	Vec3 bmin, bmax; //! Bounding box corners
	int rx, ry, rz;
	T *data;

public:
	UniformGrid(Vec3 bmin, Vec3 bmax, int res)
	: bmin(bmin), bmax(bmax) {
		rx = ry = rz = res;
		data = new T[rx*ry*rz];
	}
	UniformGrid(const UniformGrid<T>& grid) {
		bmin = grid.bmin; 
		bmax = grid.bmax;
		rx = ry = rz = grid.rx;
		data = new T[rx*ry*rz];
	}
	~UniformGrid() { delete[] data; }

	#define idx(i,j,k) ((k)*rx*ry + (j)*rx + (i))

	inline T read(int i, int j, int k) const {
		return data[idx(i,j,k)];
	}

	T interpolate(const Vec3& p) const {
		const Vec3 rel = p - bmin; // relative position
		const Vec3 q = rel.cmul(Vec3(rx-1,ry-1,rz-1)).cdiv(bmax-bmin);

		int i=(int)q.x, j=(int)q.y, k=(int)q.z;
		int i1=i+1,j1=j+1,k1=k+1;

		// If outside... return zero.
		if(i<0||j<0||k<0 || i>=rx||j>=ry||k>=rz) 
			return zero<T>();

		const Vec3 w = q - Vec3(i,j,k);
		const Vec3 w1 = Vec3(1.f,1.f,1.f) - w;
		
		return 
			read(i,j,k)    * w1.x * w1.y * w1.z +
			read(i,j,k1)   * w1.x * w1.y * w.z  +
			read(i,j1,k)   * w1.x * w.y  * w1.z +
			read(i,j1,k1)  * w1.x * w.y  * w.z  +
			read(i1,j,k)   * w.x  * w1.y * w1.z +
			read(i1,j,k1)  * w.x  * w1.y * w.z  +
			read(i1,j1,k)  * w.x  * w.y  * w1.z +
			read(i1,j1,k1) * w.x  * w.y  * w.z;
	}

	void write(int i, int j, int k, const T& val) {
		data[idx(i,j,k)] = val;
	}

	void fill(T (*func)(const Vec3&)) {
		Vec3 D = (bmax-bmin).cdiv(Vec3(rx-1,ry-1,rz-1));

		#pragma omp parallel for
		for(int i=0;i<rx;++i)
			for(int j=0;j<ry;++j)
				for(int k=0;k<rz;++k) {
					Vec3 p = Vec3(i,j,k).cmul(D) + bmin;
					data[idx(i,j,k)] = func(p);
				}
	}
	#undef idx

	Vec3 getSpacing() const {
		return (bmax-bmin).cdiv(Vec3(rx-1,ry-1,rz-1));
	}
	void getResolution(int *resX, int *resY, int *resZ) const {
		*resX = rx; 
		*resY = ry;
		*resZ = rz;
	}
};

inline void computeJacobian(const UniformGrid<Vec3>& V, UniformGrid<Mat3>& M) {
	int rx, ry, rz;
	V.getResolution(&rx,&ry,&rz);
	Vec3 H = V.getSpacing();

	#pragma omp parallel for	
	for(int i=0;i<rx;++i)
		for(int j=0;j<ry;++j)
			for(int k=0;k<rz;++k) {
				Mat3 m;

				if(i>=1 && i<rx-1) { m[0] = (V.read(i+1,j,k) - V.read(i-1,j,k)) * (1.f / (2.f * H.x) );}				
				else if(i!=0) { m[0] = (V.read(i,j,k)  - V.read(i-1,j,k)) / (H.x); }
				else { m[0] = (V.read(i+1,j,k)  - V.read(i,j,k)) / (H.x); }
				
				if(j>=1 && j<ry-1) { m[1] = (V.read(i,j+1,k) - V.read(i,j-1,k)) * (1.f / (2.f * H.y) );}				
				else if(j!=0) { m[1] = (V.read(i,j,k)  - V.read(i,j-1,k)) / (H.y); }
				else { m[1] = (V.read(i,j+1,k)  - V.read(i,j,k)) / (H.y); }

				if(k>=1 && k<rz-1) { m[2] = (V.read(i,j,k+1) - V.read(i,j,k-1)) * (1.f / (2.f * H.z) );}				
				else if(k!=0) { m[2] = (V.read(i,j,k)  - V.read(i,j,k-1)) / (H.z); }
				else { m[2] = (V.read(i,j,k+1)  - V.read(i,j,k)) / (H.z); }
				
				M.write(i,j,k,m.transpose());
			}
}

#endif
