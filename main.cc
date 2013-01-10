#define GL_GLEXT_PROTOTYPES

#include "Vec3.h"
#include "Utils.h"
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <vector>
#include "VortexSim.h"
#include "UniformGrid.h"
#include "Octree.h"
#include <GL/gl.h>

sf::RenderWindow *app;
VortexSim *sim;
std::vector<Tracer> tracers;
UniformGrid<Vec3> *velocityGrid;
UniformGrid<Mat3> *jacobianGrid;
Octree *octree;

unsigned int tracerVBO;

bool paused = false;
bool drawAxes = false;
bool drawTracers = true;
bool drawTree = false;
int width = 1920, height = 1080;

Vec3 vortonMin(-10,-10,-10), vortonMax(10,10,10);

float rand01() { return rand() * (1.f/RAND_MAX); }
Vec3 randVec3() {
	return Vec3(rand01(), rand01(), rand01()) * 2.f - Vec3(1,1,1);
}

void drawBox(const Vec3& a, const Vec3& b) {
	glBegin(GL_LINES);
	glVertex3f(a.x,a.y,a.z); glVertex3f(a.x,a.y,b.z);
	glVertex3f(a.x,a.y,a.z); glVertex3f(b.x,a.y,a.z);
	glVertex3f(a.x,a.y,b.z); glVertex3f(b.x,a.y,b.z);
	glVertex3f(b.x,a.y,a.z); glVertex3f(b.x,a.y,b.z);
	
	glVertex3f(a.x,b.y,a.z); glVertex3f(a.x,b.y,b.z);
	glVertex3f(a.x,b.y,a.z); glVertex3f(b.x,b.y,a.z);
	glVertex3f(a.x,b.y,b.z); glVertex3f(b.x,b.y,b.z);
	glVertex3f(b.x,b.y,a.z); glVertex3f(b.x,b.y,b.z);
	
	glVertex3f(a.x,a.y,a.z); glVertex3f(a.x,b.y,a.z);
	glVertex3f(b.x,a.y,b.z); glVertex3f(b.x,b.y,b.z);
	glVertex3f(b.x,a.y,a.z); glVertex3f(b.x,b.y,a.z);
	glVertex3f(a.x,a.y,b.z); glVertex3f(a.x,b.y,b.z);
	glEnd();
}

void drawAxis(const Vec3& center, const Vec3& axis) {
	Vec3 u(1,0,0), v(0,1,0);
	if(fabs(u*axis) > fabs(v*axis)) {
		v = v - (v*axis)*axis;
		v = v.normalized();
		u = axis ^ v;
		u = u.normalized();
	} else {
		u = u - (u*axis)*axis;
		u = u.normalized();
		v = axis ^ u;
		v = v.normalized();
	}	

	float radius = .01f;//axis.norm();
	int N = 6;
	glBegin(GL_LINE_LOOP);
	for(int i=0; i < N; ++i) {
		float theta = i * 3.14159265359f * 2.f / (N-1);
		float s = sinf(theta), t = cosf(theta);
		Vec3 p = center + radius * (u*s + v*t);
		glVertex3f(p.x,p.y,p.z);
	}
	glEnd();

	Vec3 axisN = axis.normalized();
	Vec3 axisA = center - axisN * .05f;
	Vec3 axisB = center + axisN * .05f;
	glBegin(GL_LINES);
	glVertex3f(axisA.x, axisA.y, axisA.z);
	glVertex3f(axisB.x, axisB.y, axisB.z);
	glEnd();
}

int frame = 0;
void glDisplay() {
	glClearColor(1,1,1,1);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30.f, (float)width/(float)height, 0.f, 100.f);
	Vec3 cam(4.1,.1,4);
	gluLookAt(cam.x,cam.y,cam.z, 0,0,0, 0,1,0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(sim==NULL) return;
	
	if(drawTree)
	{
		std::vector<Vec3> mins, maxes;
		Vec3 point = sim->particles[0].position;
		octree->getGeometryContainingPoint(point, mins, maxes);

		//octree->getGeometry(mins, maxes);
		LOG_INFO("%ld boxes collected", mins.size());
		glColor3f(0,0,0);
		for(int i=0; i<mins.size(); ++i) {
			float r = (maxes[i]-mins[i]).normSquared();
			r = .75;//1.f - expf(-1*r);
			glColor3f(r,r,r);
			drawBox(mins[i], maxes[i]);
		}
	}
	
	if(drawTracers) 
	{
		glColor4f(0,0,0, .04);
	glEnable( GL_POINT_SMOOTH );
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glPointSize( .25f );

		glBindBuffer(GL_ARRAY_BUFFER, tracerVBO);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * tracers.size(), &tracers[0], GL_STATIC_DRAW);
		glEnableClientState(GL_VERTEX_ARRAY);
		
		glVertexPointer(3, GL_FLOAT, 0, 0);
		glDrawArrays(GL_POINTS, 0, tracers.size());

		glDisableClientState(GL_VERTEX_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, 0);

		/*
	//	glPointSize(1.f);
//	glEnable( GL_POINT_SMOOTH );
//	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glPointSize( 0.5f );
		glBegin(GL_POINTS);
		for(const auto& tracer : tracers) {
			const Vec3 &p(tracer.position);

			Vec3 C(0,0,0);
			//C = tracer.density > 0 ? Vec3(1,0,0) : Vec3(0,1,0);
			glColor4f(C.x, C.y, C.z, .03f);
			glVertex3f(p.x, p.y, p.z);
		}
		glEnd();
		*/
	}

	if(drawAxes)
	{
		glLineWidth(1.f);
		glColor3f(0,.5f,0);
		for(const auto& vorton : sim->particles) {
			drawAxis(vorton.position, vorton.vorticity);
		}
	}

	glFinish();

}

inline Vec3 getVelocity(const Vec3& p) {
	//return sim->getVelocity(p);
	Vec3 f(0,0,0);
	octree->accumulateVortexForceLocalStack(p,f);	
	//octree->accumulateVortexForce(p,f);

	return f;
}

Mat3 getJacobianFromOctree(const Vec3& p) {
	float h = .04f;
	Mat3 M;
	M[0] = getVelocity(p+Vec3(h,0,0)) - getVelocity(p-Vec3(h,0,0));
	M[1] = getVelocity(p+Vec3(0,h,0)) - getVelocity(p-Vec3(0,h,0));
	M[2] = getVelocity(p+Vec3(0,0,h)) - getVelocity(p-Vec3(0,0,h));
	return M.transpose() * (1.f / (2*h));
}

void doDynamics() {
	float dt = .05f;
	sf::Clock clock;

	// Reconstruct the octree around the vortons
	clock.Reset();
	delete octree;

	vortonMin = Vec3( 1000, 1000, 1000); 
	vortonMax = Vec3(-1000,-1000,-1000);
	for(auto& vorton : sim->particles) {
		vortonMin = min(vortonMin, vorton.position);
		vortonMax = max(vortonMax, vorton.position);
	}
	vortonMin -= Vec3(.1,.1,.1);
	vortonMax += Vec3(.1,.1,.1);

	octree = new Octree( 
		(vortonMin+vortonMax)*.5f, 
		(vortonMax-vortonMin)*.5f);
	for(auto& vorton : sim->particles) 
		octree->insert(new OctreePoint(vorton.position, vorton.vorticity));
	octree->buildAggregation();
	printf("Rebuilding octree with %ld points took %ld sec\n", sim->particles.size(), clock.GetElapsedTime());

	// Compute velocity grid and its Jacobian
	//velocityGrid->fill(&getVelocity);
	//computeJacobian(*velocityGrid, *jacobianGrid);

	// Advect tracers
	clock.Reset();
	#pragma omp parallel for schedule(dynamic)
	for(int i=0; i<tracers.size(); ++i) {
		auto& tracer = tracers[i];
		tracer.position += getVelocity(tracer.position) *dt;
	}
	printf("Advected %ld tracers in %.03lf sec\n", tracers.size(), clock.GetElapsedTime());

	// Evolve vortons
	clock.Reset();
	#pragma omp parallel for
	for(int i=0; i<sim->particles.size(); ++i) {
		auto& vorton = sim->particles[i];;
		const Vec3 paddle = getVelocity(vorton.position);
		vorton.velocity = paddle;
	}
	printf("Computed velocity for %ld vortons in %.03lf sec\n", sim->particles.size(), clock.GetElapsedTime());
		
	if(1) 
	{
		clock.Reset();
		#pragma omp parallel for
		for(int i=0; i<sim->particles.size(); ++i) {
			auto& vorton = sim->particles[i];;
			const Mat3 J = getJacobianFromOctree(vorton.position);//jacobianGrid->interpolate(vorton.position);
			const Vec3 stretch = J * vorton.vorticity;
			vorton.vorticity += .25f * stretch * dt;
		}
		printf("Computed Jacobian term in %.03lf sec\n", clock.GetElapsedTime());
	}

	clock.Reset();
	#pragma omp parallel for
	for(int i=0; i<sim->particles.size(); ++i) {
		auto& vorton = sim->particles[i];
		vorton.position += vorton.velocity * dt;
	}
	printf("Advected %ld vortons in %.03lf sec\n", sim->particles.size(), clock.GetElapsedTime());

}

void takeScreenshot(char *prefix) {
	sf::Image Screen = app->Capture();
	char buff[256];
	sprintf(buff, "%s_%06d.jpg", prefix, frame);
  Screen.SaveToFile(buff);
}

void gameLoop() {
	bool running = true;
	sf::Clock clock;

	while(app->IsOpened()) {
		// Process events
		sf::Event event;
		while(app->GetEvent(event)) {
			// Window closed
			if (event.Type == sf::Event::Closed)
				app->Close();
			// Escape key pressed
			if ((event.Type == sf::Event::KeyPressed) && (event.Key.Code == sf::Key::Escape))
				app->Close();
			if ((event.Type == sf::Event::KeyPressed) && (event.Key.Code == sf::Key::Num1))
				drawTracers = !drawTracers;
			if ((event.Type == sf::Event::KeyPressed) && (event.Key.Code == sf::Key::Num2))
				drawAxes = !drawAxes;
			if ((event.Type == sf::Event::KeyPressed) && (event.Key.Code == sf::Key::Num3))
				drawTree = !drawTree;
			if ((event.Type == sf::Event::KeyPressed) && (event.Key.Code == sf::Key::Space))
				paused = !paused;
		}

		clock.Reset();

		// Dynamics
		if(!paused)
			doDynamics();

		// Display
		{
			sf::Clock displayTime;
			displayTime.Reset();
	#if 0	
			drawTree = true; drawTracers = false; drawAxes = false;
			glDisplay();
			takeScreenshot("tree");
			drawTree = false; drawTracers = false; drawAxes = true;
			glDisplay();
			takeScreenshot("vortons");
			drawTree = false; drawTracers = true; drawAxes = false;
			glDisplay();
			takeScreenshot("tracers");
			if(frame>1500) break;
	#else
			glDisplay();
	#endif
			printf("Rendered in %.03lf sec\n", displayTime.GetElapsedTime());
		}

		app->Display();
		LOG_INFO("%.02f spf / %02f fps", clock.GetElapsedTime(), 1.f/clock.GetElapsedTime());
		
		frame++;
	}
}


void init() {

	app = new sf::RenderWindow(sf::VideoMode(width, height, 32), "V o r t e x");	
	sim = new VortexSim();

	int nSegments = 256;
	const float twoPi = 2.f*M_PI;
	for(int i=0; i<nSegments; ++i) {
		float t = i * 2.f*3.14159265359f / (nSegments-1);
		Vec3 p = Vec3(cosf(t)*.5f,-2,sinf(t)*.5f) + randVec3() * .00f;

		Vec3 c = p ^ Vec3(0,1,0);
		c = c.normalized() * -3.f / nSegments;
		sim->particles.push_back(Vorton(p, c));
	}

	for(int i=0; i<nSegments; ++i) {
		float t = i * 2.f*3.14159265359f / (nSegments-1);
		Vec3 p = Vec3(cosf(t)*.45f,2,sinf(t)*.45f) + randVec3()*.00f;
		Vec3 c = p ^ Vec3(0,1,0);
		c = c.normalized() * 3.f / nSegments;
		sim->particles.push_back(Vorton(p*1.f, c));
	}

	const float dx = .01f * 2.f;
	for(float x=-.5;x<=.5;x+=dx)
	for(float y=-2;y<=2;y+=dx)
	for(float z=-.5;z<=.5;z+=dx) {
		//tracers.push_back(Tracer(Vec3(x,y+3.5,z)*.5f + randVec3()*.005f,Vec3(0,0,0)));
		//tracers.push_back(Tracer(Vec3(x,y-3.5,z)*.5f + randVec3()*.005f,Vec3(0,0,0)));
		Vec3 p = Vec3(x,y,z) * .5f;
		tracers.push_back(Tracer(p + randVec3()*.005f));
	}

	octree = new Octree(Vec3(-10,-10,-10),Vec3(10,10,10));

	// Tracer VBO
	glGenBuffers(1, &tracerVBO);
	glBindBuffer(GL_ARRAY_BUFFER, tracerVBO);

	Vec3 bmin(-4,-4,-4), bmax(4,4,4);
	int res = 100;
	velocityGrid = new UniformGrid<Vec3>(bmin, bmax, res);
	jacobianGrid = new UniformGrid<Mat3>(bmin, bmax, res);
}

int main(int argc, char **argv) {
	init();
	gameLoop();
	return 0;
}
