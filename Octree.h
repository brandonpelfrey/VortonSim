#ifndef Octree_H
#define Octree_H

#include <vector>
#include "Vec3.h"

struct OctreePoint {
  Vec3 position;
  Vec3 vorticity;
  OctreePoint(const Vec3& position)
  : position(position), vorticity(Vec3(0,0,0)) { }
  OctreePoint(const Vec3& position, const Vec3& vorticity)
  : position(position), vorticity(vorticity) { }
 
  inline Vec3 getForcing(const Vec3& p) const { 
		const float quarterPiInverse = .25f / (3.141592653589f);
		const float Vi = 1.f;

    const float radiusSquared = .01f;
		const float radiusSquaredInverse = 1.f / radiusSquared;
    const Vec3 r = p - position;
    float dist2 = .0001f + r.normSquared();
    float distInv = 1.f / sqrtf(dist2);
    float K = dist2*radiusSquaredInverse < 1.f ? distInv * radiusSquaredInverse : distInv / dist2;
    return quarterPiInverse *  (vorticity ^ r) * (Vi * K);
  }

};

class Octree {
	Vec3 origin, halfDimension;
	Vec3 bmin, bmax;
	Octree *children[8]; // Children
	OctreePoint* data; // Possible item stored at a leaf

	/*
	child:	0 1 2 3 4 5 6 7
  x:      - - - - + + + +
  y:      - - + + - - + +
	z:      - + - + - + - +
	*/

	Octree(const Octree& tree);

public:
	Octree(const Vec3& origin, const Vec3& halfDimension) 
  : origin(origin), halfDimension(halfDimension) {
		for(int i=0; i<8; ++i)
			children[i] = NULL;
    data = NULL;

		const Vec3 margin = Vec3(1,1,1) * .0f;
		bmin = origin - halfDimension - margin;
		bmax = origin + halfDimension + margin;
	}
  
  ~Octree() {
    for(int i=0; i<8; ++i) {
      if(children[i] != NULL)
        delete children[i];
    }
  }

  bool isLeafNode() const {
		/*
    for(int i=0; i<8; ++i) 
      if(children[i] != NULL)
        return false;
    return true;
		*/
		return children[0] == NULL;
  }

  int getOctantContainingPoint(const Vec3& point) const {
    int oct = 0;
    if(point.x >= origin.x) oct |= 4;
    if(point.y >= origin.y) oct |= 2;
    if(point.z >= origin.z) oct |= 1;
    return oct;
  }

  void insert(OctreePoint* point) {
    // If this node doesn't have a point yet assigned 
    // and it is a leaf, then we're done!
    if(isLeafNode()) {
      if(data==NULL) {
        data = point;
        return;
      
      } else {
        // We're at a leaf, but there's already something here
        // We will split this node and then insert the old data and this new node
  
        // Save this data point that was here for a later re-insert
        OctreePoint *oldPoint = data;

        // Split the node
		    for(int i=0; i<8; ++i) {
          // Compute new bounding box
          Vec3 newOrigin = origin;
          newOrigin.x += halfDimension.x * (i&4 ? .5f : -.5f);
          newOrigin.y += halfDimension.y * (i&2 ? .5f : -.5f);
          newOrigin.z += halfDimension.z * (i&1 ? .5f : -.5f);
          children[i] = new Octree(newOrigin, halfDimension*.5f);
        }
        
        // re-insert the old point, and insert this new point
        insert(oldPoint);
        insert(point);
      }

    } else {
      // We are at an interior node. Insert recursively into a child
      int octant = getOctantContainingPoint(point->position);
      children[octant]->insert(point);
    }
  }

  void accumulateVortexForce(const Vec3& p, Vec3& f) const {
    const Vec3 margin = Vec3(1,1,1) * 0.2 * .1;//25f;

    if(isLeafNode()) {
      if(data!=NULL) {
        f += data->getForcing(p);        
      }

    } else {
      for(int i=0; i<8; ++i) {
//        Vec3 cmin = children[i]->origin - children[i]->halfDimension - margin;
//        Vec3 cmax = children[i]->origin + children[i]->halfDimension + margin;
				const Vec3 &cmin(children[i]->bmin);
				const Vec3 &cmax(children[i]->bmax);

        bool inside = p.x<cmax.x && p.y<cmax.y && p.z<cmax.z && p.x>=cmin.x && p.y>=cmin.y && p.z>=cmin.z;

        if(inside) 
          children[i]->accumulateVortexForce(p,f);
        else if(children[i]->data != NULL) {
          f += children[i]->data->getForcing(p);
        }
      }
    }
  }

	void accumulateVortexForceLocalStack(const Vec3& p, Vec3& f) {
		Octree *todo[32];
		int stack = 0;
		const float theta = 1, thetaSquared = theta*theta;

		// push root
		todo[stack] = (this);

		while(stack>=0) {
			Octree *item = todo[stack--]; // pop item

			// Leaf?
			if(item->isLeafNode()) {
				if(item->data != NULL) {
					f += item->data->getForcing(p);
				}
			}	else {
				for(int child=0; child<8; ++child) {
					Octree *C = item->children[child];
		
					Vec3 &cmin(C->bmin), &cmax(C->bmax);
					bool inside = p.x<cmax.x && p.y<cmax.y && p.z<cmax.z && p.x>=cmin.x && p.y>=cmin.y && p.z>=cmin.z;

					if(inside && item->children[child] != NULL) {
						todo[++stack] = item->children[child];
					}
					else if(item->children[child]->data != NULL) {
						// We're not in this quadrant. See if we should recurse anway, or just
						// take the aggregate value. if distance_to_CoM * theta > L then use aggregate!
						const float L = item->children[child]->halfDimension.x * 2.f;
						float distSquared = (item->children[child]->data->position - p).normSquared();
						if(distSquared * thetaSquared > L*L) {
							f += item->children[child]->data->getForcing(p);
						}
						else {
							if(item->children[child] != NULL)
								todo[++stack] = item->children[child];
						}
					}

				}
			}
		}
	}

  void buildAggregation() {
    // Leaf nodes do nothing
    if(isLeafNode()) {
      return;
    }

    // Build children before computing things based on their results
    for(int i=0; i<8; ++i) {
      children[i]->buildAggregation();
    }

    // Compute the weighted center of mass and total vorticity
    Vec3 CoM(0,0,0), aggreateVorticity(0,0,0);
    float vorticitySum = 0.f;
    for(int i=0; i<8; ++i) {
      if(children[i]->data == NULL) continue;
      float vorticity_i = children[i]->data->vorticity.norm();
      vorticitySum += vorticity_i;
      CoM += children[i]->data->position * vorticity_i;
      aggreateVorticity += children[i]->data->vorticity;
    }

    if(vorticitySum > 1.e-5f) {
      CoM = CoM * (1.f / vorticitySum);
    }
  
    // Create aggregate  
    data = new OctreePoint(CoM, aggreateVorticity);
  }

  void getGeometry(std::vector<Vec3>& mins, std::vector<Vec3>& maxes) {
    mins.push_back(origin - halfDimension);
    maxes.push_back(origin + halfDimension);
    
    if(!isLeafNode())
    for(int i=0; i<8; ++i) {
      children[i]->getGeometry(mins,maxes);
    }
  }
  
	void getGeometryContainingPoint(const Vec3& p, std::vector<Vec3>& mins, std::vector<Vec3>& maxes) {
		if(p.x > bmax.x || p.x < bmin.x) return;
		if(p.y > bmax.y || p.y < bmin.y) return;
		if(p.z > bmax.z || p.z < bmin.z) return;

    mins.push_back(bmin);
    maxes.push_back(bmax);
    
    if(!isLeafNode())
    for(int i=0; i<8; ++i) {
      children[i]->getGeometryContainingPoint(p, mins,maxes);
    }
  }

  void getPointsInsideBox(const Vec3& bmin, const Vec3& bmax, std::vector<OctreePoint*>& results) {
    if(isLeafNode()) {
      if(data!=NULL) {
        Vec3 p = data->position;
        if(p.x>bmax.x || p.y>bmax.y || p.z>bmax.z) return;
        if(p.x<bmin.x || p.y<bmin.y || p.z<bmin.z) return;
        results.push_back(data);
      }
    } else {
      for(int i=0; i<8; ++i) {
        Vec3 cmax = children[i]->origin + children[i]->halfDimension;
        Vec3 cmin = children[i]->origin - children[i]->halfDimension;

        if(cmax.x<bmin.x || cmax.y<bmin.y || cmax.z<bmin.z) continue;
        if(cmin.x>bmax.x || cmin.y>bmax.y || cmin.z>bmax.z) continue;

        children[i]->getPointsInsideBox(bmin,bmax,results);
      } 
    }
  }

  

};

#endif
