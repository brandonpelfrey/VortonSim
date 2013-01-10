#ifndef VortexSim_H
#define VortexSim_H

#include <cmath>
#include "Vec3.h"
#include <vector>
#include <unordered_map>
#include "Utils.h"

struct Vorton {
	Vec3 position;
	Vec3 vorticity;
	Vec3 velocity;
	Vorton(const Vec3& position, const Vec3& vorticity)
	: position(position), vorticity(vorticity) { }
};

struct Tracer {
	Vec3 position;
	Tracer(const Vec3& position) 
	: position(position) { }
};

class VortexSim {

public:
	std::vector<Vorton> particles;
	VortexSim() { }

	Vec3 getVelocity(const Vec3& p) const {
		const float quarterPiInverse = .25f / (3.141592653589f);
		const float sigma = 1.e-1f;
		const float Vi = 1.f;

		Vec3 net(0,0,0);
		for(const auto &particle : particles) {
			const Vec3 r = p - particle.position;
			const float radius = r.norm();
			if(radius < sigma) {
				net += (particle.vorticity ^ r) * (Vi / (radius*sigma*sigma));
			} else {
				net += (particle.vorticity ^ r) * (Vi / (radius*radius*radius));
			}
		}
		return net * quarterPiInverse;
	}
};

#endif

