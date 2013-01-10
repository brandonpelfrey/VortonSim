#ifndef Mat3_H
#define Mat3_H

struct Mat3 {
	Vec3 rows[3];
	Mat3() { }
	Mat3(const Vec3& a, const Vec3& b, const Vec3& c) {
		rows[0] = a; rows[1] = b; rows[2] = c;
	}
	Vec3 operator*(const Vec3& v) const {
		return Vec3(rows[0]*v, rows[1]*v, rows[2]*v);
	}
	Mat3 operator+(const Mat3& m) const {
		return Mat3(rows[0]+m.rows[0], rows[1]+m.rows[1], rows[2]+m.rows[2]);
	}
	Mat3 operator*(const float r) const {
		return Mat3(rows[0]*r, rows[1]*r, rows[2]*r);
	}

	Mat3 transpose() const {
		Mat3 R;
		R.rows[0] = Vec3(rows[0].x, rows[1].x, rows[2].x);
		R.rows[1] = Vec3(rows[0].y, rows[1].y, rows[2].y);
		R.rows[2] = Vec3(rows[0].z, rows[1].z, rows[2].z);
		return R;
	}

	Vec3& operator[](const unsigned int i) { return rows[i]; }
	const Vec3& operator[](const unsigned int i) const { return rows[i]; }
};

#endif
