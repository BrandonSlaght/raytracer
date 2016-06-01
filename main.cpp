#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
 
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#define NOMINMAX
#include <windows.h>
#endif // Win32 platform
 
#include <GL/gl.h>
#include <GL/glu.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/
#include <GL/glut.h>

#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>

// simple material class, with object color, and headlight shading
class Material
{
protected:
	float3 color;
public:
	
	Material(float3 color) :color(color) {}

	virtual float3 getReflectionDir(float3 inDir, float3 normal) {
		return float3(0, 0, 0);
	}

	virtual float3 getRefractionDir(float3 inDir, float3 normal) {
		return float3(0, 0, 0);
	}

	virtual float3 getReflectance() {
		return float3(0, 0, 0);
	}

	virtual float getRefractance() {
		return 0;
	}

	virtual float3 shade(
		float3 position,
		float3 normal,
		float3 viewDir,
		float3 lightDir,
		float3 lightPowerDensity)
	{
		return float3(0, 0, 0);
	}

	virtual float3 getColor(
		float3 position,
		float3 normal,
		float3 viewDir)
	{
		float f = viewDir.dot(normal);
		return f >= 0 ? color*f : float3(0.0, 0.0, 0.0);
	}

	float snoise(float3 r) {
		unsigned int x = 0x0625DF73;
		unsigned int y = 0xD1B84B45;
		unsigned int z = 0x152AD8D0;
		float f = 0;
		for (int i = 0; i<32; i++) {
			float3 s(x / (float)0xffffffff,
				y / (float)0xffffffff,
				z / (float)0xffffffff);
			f += sin(s.dot(r));
			x = x << 1 | x >> 31;
			y = y << 1 | y >> 31;
			z = z << 1 | z >> 31;
		}
		return f / 64.0 + 0.5;
	}
};

class DiffuseMaterial : public Material
{
public:
	DiffuseMaterial(float3 m) : 
		Material(m)
		{}

	float3 shade(
		float3 position,
		float3 normal,
		float3 viewDir,
		float3 lightDir,
		float3 lightPowerDensity)
	{
		float cosTheta = normal.dot(lightDir);
		if (cosTheta < 0) return float3(0, 0, 0);
		return color * lightPowerDensity * cosTheta;
	}
};

class SolidMaterial : public DiffuseMaterial {
	float3 frontFaceColor;
	float3 backFaceColor;
public:

	SolidMaterial(float3 frontFaceColor, float3 backFaceColor) :
		DiffuseMaterial(frontFaceColor),
		frontFaceColor(frontFaceColor),
		backFaceColor(backFaceColor) {}

	float3 shade(
		float3 position,
		float3 normal,
		float3 viewDir,
		float3 lightDir,
		float3 lightPowerDensity)
	{
		float cosTheta = normal.dot(lightDir);
		if (cosTheta < 0) return float3(0, 0, 0);
		float f = viewDir.dot(normal);
		if (f >= 0)
			return frontFaceColor * lightPowerDensity * cosTheta;
		else
			return backFaceColor * lightPowerDensity * cosTheta;
	}

	//virtual float3 getColor(
	//	float3 position, float3 normal, float3 viewDir)
	//{
	//	float f = viewDir.dot(normal);
	//	if (f >= 0)
	//		return frontFaceColor*f;
	//	else
	//		return backFaceColor*f;
	//}
};

class Wood : public DiffuseMaterial
{
	float scale;
	float turbulence;
	float period;
	float sharpness;
public:
	Wood() :
		DiffuseMaterial(float3(1, 1, 1))
	{
		scale = 16;
		turbulence = 500;
		period = 8;
		sharpness = 10;
	}
	virtual float3 shade(
		float3 position,
		float3 normal,
		float3 viewDir,
		float3 lightDir,
		float3 lightPowerDensity)
	{
		//return normal;
		float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence + 10000.0;
		w -= int(w);
		return ((float3(1, 0.3, 0) * w + float3(0.35, 0.1, 0.05) * (1 - w)) * normal.dot(viewDir));
	}
};

class PhongMaterial : public DiffuseMaterial {
	float shininess;
public:
	PhongMaterial(float3 material, float shininess) :
		DiffuseMaterial(material),
		shininess(shininess) 
	{
	}

	float3 shade(float3 position, float3 normal, float3 viewDir, float3 lightDir, float3 lightPowerDensity)
	{
		float cosTheta = normal.dot(lightDir);
		if (cosTheta < 0) return float3(0, 0, 0);
		float3 halfway = (viewDir + lightDir).normalize();
		float cosDelta = normal.dot(halfway);
		if (cosDelta < 0) return float3(0, 0, 0);
		return lightPowerDensity * color * pow(cosDelta, shininess);
	}

};

class IdealAbsorbingMaterial : public DiffuseMaterial {
	float3  reflectance;
public:

	IdealAbsorbingMaterial(float3 color, float3 reflectance) : DiffuseMaterial(color), reflectance(reflectance) {}

	float3 getReflectionDir(float3 inDir, float3 normal) {
		float3 perp = -normal * normal.dot(inDir);
		float3 parallel = inDir - perp;
		return parallel - perp;
	}

	float3 getReflectance() {
		return reflectance;
	}
};

class Metal : public IdealAbsorbingMaterial
{
	float scale;
	float turbulence;
	float period;
	float sharpness;
public:
	Metal() :
		IdealAbsorbingMaterial(float3(.7, .7, .7), float3(.5, .5, .5))
	{
		scale = 32;
		turbulence = 50;
		period = 32;
		sharpness = 1;
	}
	virtual float3 shade(
		float3 position,
		float3 normal,
		float3 viewDir,
		float3 lightDir,
		float3 lightPowerDensity)
	{
		//return normal;
		float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence;
		w = pow(sin(w)*0.5 + 0.5, 4);
		return (float3(.5, .5, .5) * w + float3(.6, .6, .6) * (1 - w)) * normal.dot(viewDir);
	}
};

class IdealTransparentMaterial : public DiffuseMaterial {
	float  refractiveIndex;
	float  refractance;
public:
	IdealTransparentMaterial(float3 color, float refractiveIndex, float  refractance) : DiffuseMaterial(color), refractiveIndex(refractiveIndex), refractance(refractance) {}

	float3 getRefractionDir(float3 inDir, float3 normal) {
		float ri = refractiveIndex;
		float cosa = -normal.dot(inDir);
		if (cosa < 0) { 
			cosa = -cosa; 
			normal = -normal; 
			ri = 1 / ri; 
		}
		float disc = 1 - (1 - cosa * cosa) / ri / ri;
		float cosb = (disc < 0) ? 0 : sqrt(disc);
		float3 perp = -normal * normal.dot(inDir);
		float3 parallel = inDir - perp;
		//return parallel;
		return parallel * (1.0/ri) - normal * cosb;
	}

	float getRefractance() {
		return refractance;
	}
};

class LightSource
{
public:
	virtual float3 getPowerDensityAt(float3 x) = 0;
	virtual float3 getLightDirAt(float3 x) = 0;
	virtual float  getDistanceFrom(float3 x) = 0;
};

class DirectionalLight : public LightSource
{
	float power;
	float3 direction;
public:
	DirectionalLight(float power, float3 direction) :
		power(power),
		direction(direction) 
	{
	}
	float3 getPowerDensityAt(float3 x) {
		return float3(power, power, power);
	}

	float3 getLightDirAt(float3 x) {
		return direction;
	}

	float getDistanceFrom(float3 x) {
		return FLT_MAX;
	}
};

class PointLight : public LightSource 
{
	float3 power;
	float3 position;
public:
	PointLight(float3 power, float3 position) :
		power(power),
		position(position) 
	{
	}
	float3 getPowerDensityAt(float3 x) {
		return power*(1.0 / (4.0 * M_PI*((position - x).dot(position - x)))); //vector.itself = vector size ^2  
	}
	float3 getLightDirAt(float3 x) {
		return (position - x).normalize();
	}
	float getDistanceFrom(float3 x) {
		return sqrtf(x.x*position.x + x.y*position.y + x.z*position.z);
	}
};

// Skeletal camera class.
class Camera
{
	float3 eye;		//< world space camera position
	float3 lookAt;	//< center of window in world space
	float3 right;	//< vector from window center to window right-mid (in world space)
	float3 up;		//< vector from window center to window top-mid (in world space)

public:
	Camera()
	{
		eye = float3(0, 0, 3);
		lookAt = float3(0, 0, 2);
		right = float3(1, 0, 0);
		up = float3(0, 1, 0);
	}
	float3 getEye()
	{
		return eye;
	}
	// compute ray through pixel at normalized device coordinates
	float3 rayDirFromNdc(const float2 ndc) {
		return (lookAt - eye
			+ right * ndc.x
			+ up    * ndc.y
			).normalize();
	}
};

// Ray structure.
class Ray
{
public:
    float3 origin;
    float3 dir;
    Ray(float3 o, float3 d)
    {
        origin = o;
        dir = d;
    }
};

// Hit record structure. Contains all data that describes a ray-object intersection point.
class Hit
{
public:
	Hit()
	{
		t = -1;
	}
	float t;				//< Ray paramter at intersection. Negative means no valid intersection.
	float3 position;		//< Intersection coordinates.
	float3 normal;			//< Surface normal at intersection.
	Material* material;		//< Material of intersected surface.
};

// Object abstract base class.
class Intersectable
{
protected:
	Material* material;
public:
	Intersectable(Material* material):material(material) {}
    virtual Hit intersect(const Ray& ray)=0;
};

// Simple helper class to solve quadratic equations with the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and store the results.
class QuadraticRoots
{
public:
	float t1;
	float t2;
	// Solves the quadratic a*t*t + b*t + c = 0 using the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and set members t1 and t2 to store the roots.
	QuadraticRoots(float a, float b, float c)
	{
        float discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) // no roots
		{
			t1 = -1;
			t2 = -1;
			return;
		}
        float sqrt_discr = sqrt( discr );
		t1 = (-b + sqrt_discr)/2.0/a;
		t2 = (-b - sqrt_discr)/2.0/a;
	}
	// Returns the lesser of the positive solutions, or a negative value if there was no positive solution.
	float getLesserPositive()
	{
		return (0 < t1 && t1 < t2 || t2 < 0)?t1:t2;
	}
};

// Object realization.
class Sphere : public Intersectable
{
	float3 center;
	float radius;
public:
    Sphere(const float3& center, float radius, Material* material):
		Intersectable(material),
		center(center),
		radius(radius)
    {
    }
	QuadraticRoots solveQuadratic(const Ray& ray)
	{
        float3 diff = ray.origin - center;
        float a = ray.dir.dot(ray.dir);
        float b = diff.dot(ray.dir) * 2.0;
        float c = diff.dot(diff) - radius * radius;
		return QuadraticRoots(a, b, c);
	}
	float3 getNormalAt(float3 r)
	{
		return (r - center).normalize();
	}
    Hit intersect(const Ray& ray)
    {
		// This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal
		float t = solveQuadratic(ray).getLesserPositive();
			
		Hit hit;
		hit.t = t;
		hit.material = material;
		hit.position = ray.origin + ray.dir * t;
		hit.normal = getNormalAt(hit.position);

		return hit;
    }
}; 

// CLASS PLANE COULD COME HERE
class Plane : public Intersectable
{
	float3 center;
	float3 normal;
public:
	Plane(const float3& center, float3& normal, Material* material) :
		Intersectable(material),
		center(center),
		normal(normal)
	{
	}
	float3 getNormalAt(float3 r)
	{
		return normal;
	}
	Hit intersect(const Ray& ray)
	{
		// This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal
		float t = normal.dot(center - ray.origin) / normal.dot(ray.dir);

		Hit hit;
		hit.t = t;
		hit.material = material;
		hit.position = ray.origin + ray.dir * t;
		hit.normal = getNormalAt(hit.position);

		return hit;
	}
};

// CLASS QUADRIC COULD COME HERE
class Quadric : public Intersectable
{
	float4x4 coeffs;
public:
	Quadric(float4x4 coeffs, Material* material) :
		Intersectable(material),
		coeffs(coeffs)
	{
	}
	QuadraticRoots solveQuadratic(const Ray& ray)
	{
		float4 direction = float4(ray.dir);
		direction.w = 0.0;
		float4 origin = float4(ray.origin);
		float a = direction.dot(coeffs * direction);
		float b = direction.dot(coeffs * origin) + origin.dot(coeffs * direction);
		float c = origin.dot(coeffs * origin);
		return QuadraticRoots(a, b, c);
	}
	float3 getNormalAt(float3 r)
	{
		float4 r4 = float4(r);
		float4 calc = coeffs*r4 + r4*coeffs;
		return calc.drop().normalize();
	}
	Hit intersect(const Ray& ray)
	{
		// This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal
		float t = solveQuadratic(ray).getLesserPositive();

		Hit hit;
		hit.t = t;
		hit.material = material;
		hit.position = ray.origin + ray.dir * t;
		hit.normal = getNormalAt(hit.position);

		return hit;
	}
	void transform(float4x4 transformation) {
		coeffs = transformation.invert()*coeffs*transformation.invert().transpose();
	}
	Quadric* parallelPlanes() {
		float4x4 A = float4x4::identity();
		A._00 = 0;
		A._11 = 1;
		A._22 = 0;
		A._33 = -1;
		return this;
	}
	bool contains(float3 r) {
		float rhomo = float4(r).dot(coeffs * float4(r));
		return rhomo > 0 ? false : true;
	}
};

// CLASS CLIPPEDQUADRIC COULD COME HERE
class Clipped : public Intersectable
{
	Quadric shape;
	std::vector<Quadric*> clippers;
	std::vector<Quadric*> anticlippers;
public:
	Clipped(Quadric shape, Material* material) :
		Intersectable(material),
		shape(shape)
	{
	}
	void addClipper(Quadric* clipper) {
		clippers.push_back(clipper);
	}
	void addAnticlipper(Quadric* clipper) {
		anticlippers.push_back(clipper);
	}
	QuadraticRoots solveQuadratic(const Ray& ray)
	{
		return shape.solveQuadratic(ray);
	}
	Hit intersect(const Ray& ray)
	{			
		float st1 = solveQuadratic(ray).t1;
		float st2 = solveQuadratic(ray).t2;
		Hit hit1;
		hit1.t = st1;
		hit1.material = material;
		hit1.position = ray.origin + ray.dir * st1;
		hit1.normal = shape.getNormalAt(hit1.position);
		Hit hit2;
		hit2.t = st2;
		hit2.material = material;
		hit2.position = ray.origin + ray.dir * st2;
		hit2.normal = shape.getNormalAt(hit2.position);
		for (int i = 0; i < clippers.size(); i++) {
			!clippers[i]->contains(hit1.position) ? hit1.t = -1 : true;
			!clippers[i]->contains(hit2.position) ? hit2.t = -1 : true;	
		}
		for (int i = 0; i < anticlippers.size(); i++) {
			anticlippers[i]->contains(hit1.position) ? hit1.t = -1 : true;
			anticlippers[i]->contains(hit2.position) ? hit2.t = -1 : true;
		}
		if ((0 < hit1.t && hit1.t < hit2.t) || hit2.t < 0) {
			return hit1;
		}
		//if ((0 < hit2.t && hit2.t < hit1.t) || hit1.t < 0) {
		//	return hit2;
		//}
		return hit2;

		//return (0 < hit1.t && hit1.t < hit2.t || hit2.t < 0) ? hit1 : hit2;
	}
	void transform(float4x4 transformation) {
		shape.transform(transformation);
		for (int i = 0; i < clippers.size(); i++) {
			clippers[i]->transform(transformation);
		}
		for (int i = 0; i < anticlippers.size(); i++) {
			anticlippers[i]->transform(transformation);
		}
	}
};

class Scene
{
	Camera camera;
	std::vector<Intersectable*> objects;
	std::vector<Material*> materials;

	std::vector<LightSource*> lights;
public:
	Scene()
	{
		//Material* m1 = new Material(float3(1, 1, 0));
		//Material* m2 = new Material(float3(0, 1, 1));
		Material* m3 = new IdealTransparentMaterial(float3(1, 0, 1), .5, .5);
		//Material* m4 = new IdealAbsorbingMaterial(float3(1, 0, 0), float3(1, 1, 1));
		//Material* m5 = new SolidMaterial(float3(1, 1, 0), float3(0, 1, 1));
		//Material* m6 = new SolidMaterial(float3(1, 1, 0), float3(1, 1, 0));
		//Material* m7 = new Wood();
		//Material* m8 = new Marble();
		//Material* m9 = new DiffuseMaterial(float3(1, 0, 1));
		//Material* m10 = new DiffuseMaterial(float3(1, 1, 0));
		//Material* m11 = new PhongMaterial(float3(1, 1, 0), .5);
		//Material* m12 = new PhongMaterial(float3(1, 1, 0), .1);
		//materials.push_back(m1);
		//materials.push_back(m2);
		//materials.push_back(m3);
		//materials.push_back(m4);
		//materials.push_back(m5);
		//materials.push_back(m6);
		//materials.push_back(m7);
		//materials.push_back(m8);
		//materials.push_back(m9);
		//materials.push_back(m10);
		//materials.push_back(m11);
		//materials.push_back(m12);
		//Sphere* s1 = new Sphere(float3(0.5, 0.5, 0.1), .5, m1); //materials
		//objects.push_back(s1);
		//Sphere* s2 = new Sphere(float3(0.5, -0.5, 0.1), .5, m2);
		//objects.push_back(s2);
		Sphere* s3 = new Sphere(float3(0, -.5, 0), 1, m3); //idealabsorbingmaterials
		//objects.push_back(s3);
		//Sphere* s4 = new Sphere(float3(-0.5, -0.5, 0.1), .5, m4);
		//objects.push_back(s4);
		//Sphere* s5 = new Sphere(float3(1.5, .5, 0.1), .5, m5); //solidmaterials
		//objects.push_back(s5);
		//Sphere* s6 = new Sphere(float3(1.5, -0.5, 0.1), .5, m6);
		//objects.push_back(s6);
		//Sphere* s7 = new Sphere(float3(-1.5, -0.5, 0.1), .5, m7); //noise materials
		//objects.push_back(s7);
		//Sphere* s8 = new Sphere(float3(-1.5, 0.5, 0.1), .5, m8);
		//objects.push_back(s8);
		//Sphere* s9 = new Sphere(float3(-1.5, 1.5, 0.1), .5, m9); //diffusematerials
		//objects.push_back(s9);
		//Sphere* s10 = new Sphere(float3(-0.5, 1.5, 0.1), .5, m10);
		//objects.push_back(s10);
		//Sphere* s11 = new Sphere(float3(0.5, 1.5, 0.1), .5, m11); //phongmaterials 
		//objects.push_back(s11);
		//Sphere* s12 = new Sphere(float3(1.5, 1.5, 0.1), .5, m12);
		//objects.push_back(s12);
		//Plane* p1 = new Plane(float3(0.2, -1, 0.1), float3(0, 1, 0.1), m3);
		//objects.push_back(p1);
		//float4x4 no = float4x4().identity();
		//no.m[3][3] = -.1;
		//Quadric q1 = Quadric(no, m9);
		//Quadric q2 = Quadric(no, m10);
		//Quadric q3 = Quadric(no, m9);
		//q3.transform(float4x4::scaling(float3(1, 1.5, 1)));
		//q2.transform(float4x4::scaling(float3(1.5, 1, 1)));
		////objects.push_back(q3);
		////q2.transform(no.translation(float3(1, -1, 0)));
		//Clipped* c1 = new Clipped(q2, q3, m10);
		////objects.push_back(c1);
		////Quadric q4 = Quadric(q3->transform(float4x4::scaling(float3(1, 2, 1))* float4x4::rotation(float3(0, 0, 1), 1)),m10);

		//LIGHTS
		DirectionalLight* left = new DirectionalLight(.1, float3(-2, 1, 3));
		lights.push_back(left);
		DirectionalLight* right = new DirectionalLight(.1, float3(2, 1, 3));
		lights.push_back(right);
		DirectionalLight* d3 = new DirectionalLight(.5, float3(0, 0, 3));
		//lights.push_back(d3);
		PointLight* pl1 = new PointLight(float3(50, 50, 50), float3(0, 0, 3));
		//lights.push_back(pl1);

		//STAGE
		Material* wood = new Wood();
		materials.push_back(wood);
		float4x4 p1 = float4x4().identity();
		p1.m[0][0] = 0;
		p1.m[2][2] = 0;
		p1.m[3][3] = -1;
		float4x4 p2 = float4x4().identity();
		p2.m[0][0] = 0;
		p2.m[2][2] = 0;
		p2.m[3][3] = -1;
		Plane* p0 = new Plane(float3(0.2, -1, 0.1), float3(0, 1, 0.1), wood);
		Quadric q = Quadric(p1, wood);
		Quadric* q0 = new Quadric(p2, wood);
		q0->transform(float4x4::rotation(float3(1, .5, 1), 90));
		Clipped* c = new Clipped(q, wood);
		c->addClipper(q0);
		objects.push_back(p0);

		//DRUM
		Material* drum = new IdealAbsorbingMaterial(float3(.9, .9, .9), float3(.9, .9, .9));
		float4x4 drumcylinder = float4x4().identity();
		drumcylinder.m[1][1] = 0;
		drumcylinder.m[3][3] = -.1;
		float4x4 drumclip = float4x4().identity();
		drumclip.m[0][0] = 0;
		drumclip.m[2][2] = 0;
		drumclip.m[3][3] = -1;
		Quadric drumq1 = Quadric(drumcylinder, drum);
		Quadric* drumq2 = new Quadric(drumclip, drum);
		Clipped* drumobject = new Clipped(drumq1, drum);
		drumobject->addClipper(drumq2);
		drumobject->transform(float4x4::scaling(float3(1.2, .5, 1.2)));
		drumobject->transform(float4x4::translation(float3(-2, -.5, 0)));
		objects.push_back(drumobject);

		//FLUTE
		Material* flute = new DiffuseMaterial(float3(.8, .8, .8));
		float4x4 flutecylinder = float4x4().identity();
		flutecylinder.m[1][1] = 0;
		flutecylinder.m[3][3] = -.1;
		float4x4 fluteclipper = float4x4().identity();
		fluteclipper.m[0][0] = 0;
		fluteclipper.m[2][2] = 0;
		fluteclipper.m[3][3] = -1;
		Quadric fluteq1 = Quadric(flutecylinder, flute);
		Quadric* fluteq2 = new Quadric(fluteclipper, flute);
		fluteq2->transform(float4x4::scaling(float3(1.5, 1.5, 1.5)));
		Clipped* fluteobject = new Clipped(fluteq1, flute);
		fluteobject->addClipper(fluteq2);

		float4x4 flutehole1 = float4x4().identity();
		flutehole1.m[1][1] = 0;
		flutehole1.m[3][3] = -.1;
		Quadric* fluteholeq1 = new Quadric(flutehole1, flute);
		fluteholeq1->transform(float4x4::rotation(float3(0, 0, 1), 1.5));
		fluteholeq1->transform(float4x4::scaling(float3(.5, .5, .5)));
		fluteholeq1->transform(float4x4::translation(float3(0, .75, 0)));
		fluteobject->addAnticlipper(fluteholeq1);

		float4x4 flutehole2 = float4x4().identity();
		flutehole2.m[1][1] = 0;
		flutehole2.m[3][3] = -.1;
		Quadric* fluteholeq2 = new Quadric(flutehole2, flute);
		fluteholeq2->transform(float4x4::rotation(float3(0, 0, 1), 1.5));
		fluteholeq2->transform(float4x4::scaling(float3(.5, .5, .5)));
		fluteholeq2->transform(float4x4::translation(float3(0, .25, 0)));
		fluteobject->addAnticlipper(fluteholeq2);

		float4x4 flutehole3 = float4x4().identity();
		flutehole3.m[1][1] = 0;
		flutehole3.m[3][3] = -.1;
		Quadric* fluteholeq3 = new Quadric(flutehole3, flute);
		fluteholeq3->transform(float4x4::rotation(float3(0, 0, 1), 1.5));
		fluteholeq3->transform(float4x4::scaling(float3(.5, .5, .5)));
		fluteholeq3->transform(float4x4::translation(float3(0, -.25, 0)));
		fluteobject->addAnticlipper(fluteholeq3);

		float4x4 flutehole4 = float4x4().identity();
		flutehole4.m[1][1] = 0;
		flutehole4.m[3][3] = -.1;
		Quadric* fluteholeq4 = new Quadric(flutehole4, flute);
		fluteholeq4->transform(float4x4::rotation(float3(0, 0, 1), 1.5));
		fluteholeq4->transform(float4x4::scaling(float3(.5, .5, .5)));
		fluteholeq4->transform(float4x4::translation(float3(0, -.75, 0)));
		fluteobject->addAnticlipper(fluteholeq4);

		fluteobject->transform(float4x4::scaling(float3(.5, 1, .5)));
		fluteobject->transform(float4x4::rotation(float3(.2, .4, .5), 55));
		fluteobject->transform(float4x4::translation(float3(-2, 1.3, -2)));
		objects.push_back(fluteobject);

		//TRUMPET
		Material* trumpet = new DiffuseMaterial(float3(.8, .8, .1));
		float4x4 trumpetcylinder = float4x4().identity();
		trumpetcylinder.m[1][1] = 0;
		trumpetcylinder.m[3][3] = -.1;
		float4x4 trumpetclipper1 = float4x4().identity();
		trumpetclipper1.m[0][0] = 0;
		trumpetclipper1.m[2][2] = 0;
		trumpetclipper1.m[3][3] = -1;
		Quadric trumpetq1 = Quadric(trumpetcylinder, trumpet);
		Quadric* trumpetq2 = new Quadric(trumpetclipper1, trumpet);
		Clipped* trumpetobject1 = new Clipped(trumpetq1, trumpet);
		trumpetobject1->addClipper(trumpetq2);
		trumpetobject1->transform(float4x4::scaling(float3(.4, 1, .4)));
		trumpetobject1->transform(float4x4::rotation(float3(.2, .4, -.5), 55));
		trumpetobject1->transform(float4x4::translation(float3(-1.5, 0, -1)));
		objects.push_back(trumpetobject1);

		float4x4 hyperboloid = float4x4().identity();
		hyperboloid.m[0][0] = 1;
		hyperboloid.m[1][1] = -1;
		hyperboloid.m[2][2] = 1;
		hyperboloid.m[3][3] = -.1;
		float4x4 trumpetclipper2 = float4x4().identity();
		trumpetclipper2.m[0][0] = 0;
		trumpetclipper2.m[2][2] = 0;
		trumpetclipper2.m[3][3] = -1;
		Quadric trumpetq3 = Quadric(hyperboloid, trumpet);
		Quadric* trumpetq4 = new Quadric(trumpetclipper2, trumpet);
		trumpetq4->transform(float4x4::scaling(float3(.5, .5, .5)));
		trumpetq4->transform(float4x4::translation(float3(0, -.5, 0)));
		Clipped* trumpetobject2 = new Clipped(trumpetq3, trumpet);
		trumpetobject2->addClipper(trumpetq4);
		trumpetobject2->transform(float4x4::translation(float3(0, -.5, 0)));
		trumpetobject2->transform(float4x4::scaling(float3(.4, 1, .4)));
		trumpetobject2->transform(float4x4::rotation(float3(.2, .4, -.5), 55));
		trumpetobject2->transform(float4x4::translation(float3(-1.5, 0, -1)));
		objects.push_back(trumpetobject2);

		//MARACCAS
		Material* maraca = new PhongMaterial(float3(0, .7, 1), -1);
		float4x4 cylinder3 = float4x4().identity();
		cylinder3.m[1][1] = 0;
		cylinder3.m[3][3] = -.1;
		float4x4 plane3 = float4x4().identity();
		plane3.m[0][0] = 0;
		plane3.m[2][2] = 0;
		plane3.m[3][3] = -1;
		Quadric q7 = Quadric(cylinder3, maraca);
		Quadric* q8 = new Quadric(plane3, maraca);
		Clipped* c4 = new Clipped(q7, maraca);
		c4->addClipper(q8);
		c4->transform(float4x4::scaling(float3(.2, 1, .2)));
		c4->transform(float4x4::rotation(float3(.2, 1.1, .8), 55));
		c4->transform(float4x4::translation(float3(1, 0, -1)));
		objects.push_back(c4);
		
		float4x4 oval = float4x4().identity();
		oval.m[3][3] = -.1;
		Quadric* q9 = new Quadric(oval, maraca);
		q9->transform(float4x4::translation(float3(0, .5, 0)));
		q9->transform(float4x4::scaling(float3(1, 1.8, 1)));
		q9->transform(float4x4::rotation(float3(.2, 1.1, .8), 55));
		q9->transform(float4x4::translation(float3(1, 0, -1)));
		objects.push_back(q9);

		float4x4 cylinder4 = float4x4().identity();
		cylinder4.m[1][1] = 0;
		cylinder4.m[3][3] = -.1;
		float4x4 plane4 = float4x4().identity();
		plane4.m[0][0] = 0;
		plane4.m[2][2] = 0;
		plane4.m[3][3] = -1;
		Quadric q10 = Quadric(cylinder4, maraca);
		Quadric* q11 = new Quadric(plane4, maraca);
		Clipped* c5 = new Clipped(q10, maraca);
		c5->addClipper(q11);
		c5->transform(float4x4::scaling(float3(.2, 1, .2)));
		c5->transform(float4x4::rotation(float3(.2, -1.1, -.8), 55));
		c5->transform(float4x4::translation(float3(1, 0, -1.2)));
		objects.push_back(c5);
	
		Quadric* q12 = new Quadric(oval, maraca);
		q12->transform(float4x4::translation(float3(0, .5, 0)));
		q12->transform(float4x4::scaling(float3(1, 1.8, 1)));
		q12->transform(float4x4::rotation(float3(.2, -1.1, -.8), 55));
		q12->transform(float4x4::translation(float3(1, 0, -1.2)));
		objects.push_back(q12);

		//BELL
		Material* bell = new Metal();
		float4x4 paraboloid = float4x4().identity();
		paraboloid.m[0][0] = 1;
		paraboloid.m[1][1] = -2;
		paraboloid.m[2][2] = 1;
		paraboloid.m[3][3] = 1;
		float4x4 bellclipper = float4x4().identity();
		bellclipper.m[0][0] = 0;
		bellclipper.m[2][2] = 0;
		bellclipper.m[3][3] = -1;
		Quadric bellq1 = Quadric(paraboloid, bell);
		Quadric* bellq2 = new Quadric(bellclipper, bell);
		bellq2->transform(float4x4::translation(float3(0,-1,0)));
		Clipped* bellobject = new Clipped(bellq1, bell);
		bellobject->addClipper(bellq2);
		bellobject->transform(float4x4::scaling(float3(.2, .9, .2)));
		bellobject->transform(float4x4::rotation(float3(.2, 1.1, .8), 5));
		bellobject->transform(float4x4::translation(float3(-1, 3.5, -1.5)));
		objects.push_back(bellobject);

		//VIOLIN
		Material* violin = new DiffuseMaterial(float3(.6, .3, .3));
		float4x4 neck = float4x4().identity();
		neck.m[1][1] = 0;
		neck.m[3][3] = -.1;
		float4x4 neckclipper = float4x4().identity();
		neckclipper.m[0][0] = 0;
		neckclipper.m[2][2] = 0;
		neckclipper.m[3][3] = -1;
		Quadric neckq1 = Quadric(neck, violin);
		Quadric* neckq2 = new Quadric(neckclipper, violin);
		Clipped* neckobject = new Clipped(neckq1, violin);
		neckobject->addClipper(neckq2);
		neckobject->transform(float4x4::scaling(float3(.3, 1, .1)));
		neckobject->transform(float4x4::scaling(float3(.75, .75, .75)));
		neckobject->transform(float4x4::translation(float3(2.5, 1.5, -1)));
		neckobject->transform(float4x4::rotation(float3(.1, -.5, 0), .2));
		objects.push_back(neckobject);

		float4x4 body = float4x4().identity();
		body.m[3][3] = -.1;
		float4x4 bodyclipper1 = float4x4().identity();
		bodyclipper1.m[3][3] = -.1;
		Quadric bodyq1 = Quadric(body, violin);
		bodyq1.transform(float4x4::scaling(float3(2, 4, .3)));
		Quadric* bodyq2 = new Quadric(bodyclipper1, violin);
		bodyq2->transform(float4x4::translation(float3(0, -.35, .3)));

		float4x4 bodyclipper2 = float4x4().identity();
		bodyclipper2.m[3][3] = -.1;
		Quadric* bodyq3 = new Quadric(bodyclipper2, violin);
		bodyq3->transform(float4x4::scaling(float3(2, 2, 2)));
		bodyq3->transform(float4x4::translation(float3(1, .25, 0)));

		float4x4 bodyclipper3 = float4x4().identity();
		bodyclipper3.m[3][3] = -.1;
		Quadric* bodyq4 = new Quadric(bodyclipper3, violin);
		bodyq4->transform(float4x4::scaling(float3(2, 2, 2)));
		bodyq4->transform(float4x4::translation(float3(-1, .25, 0)));

		Clipped* bodyobject = new Clipped(bodyq1, violin);
		bodyobject->addAnticlipper(bodyq2);
		bodyobject->addAnticlipper(bodyq3);
		bodyobject->addAnticlipper(bodyq4);
		bodyobject->transform(float4x4::scaling(float3(.75, .75, .75)));
		bodyobject->transform(float4x4::translation(float3(2.5, 0, -1)));
		bodyobject->transform(float4x4::rotation(float3(.1, -.5, 0), .2));
		objects.push_back(bodyobject);
	}

	~Scene()
	{
		// UNCOMMENT THESE WHEN APPROPRIATE
		for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
			delete *iMaterial;
		for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
			delete *iObject;		
		for (std::vector<LightSource*>::iterator iLight = lights.begin(); iLight != lights.end(); ++iLight)
			delete *iLight;
	}

public:
	Camera& getCamera()
	{
		return camera;
	}
	
	Hit getFirstIntersect(Ray ray)
	{
		Hit h1;
		h1.t = FLT_MAX;

		for (int i = 0; i < objects.size(); i++) {
			//printf("%i\n", i);
			Hit h2 = objects[i]->intersect(ray);
			if (h2.t < h1.t && h2.t > 0) {
				h1 = h2;
			}
		}
		return h1;
	}

	float3 trace(const Ray& ray, int depth)
	{

		Hit h1 = getFirstIntersect(ray);

		if(h1.t <= 0 || h1.t == FLT_MAX)
			return float3(0, 0, 0);

		float3 shade = float3(0, 0, 0);

		for (int i = 0; i < lights.size(); i++) {
			Ray shadowRay = Ray(h1.position + h1.normal*.01, lights[i]->getLightDirAt(h1.position));
			Hit shadow = getFirstIntersect(shadowRay);
			if (shadow.t >= lights[i]->getDistanceFrom(h1.position)) {
				shade += h1.material->shade(h1.position, h1.normal, -ray.dir, lights[i]->getLightDirAt(h1.position), lights[i]->getPowerDensityAt(h1.position));
			}
		}		

		float3 reflect = h1.material->getReflectionDir(-ray.dir, h1.normal);
		Ray reflectRay = Ray(h1.position + h1.normal*.01, reflect);
		if (depth > 0) {
			float before = shade.x;
			shade += h1.material->getReflectance() * trace(reflectRay, depth - 1);
			if (before != shade.x) {
				//printf("b4 : %f\n\n", before);
				//printf("af : %f\n\n", shade.x);
			}
		}

		float3 refract = h1.material->getRefractionDir(-ray.dir, h1.normal);
		Ray refractRay = Ray(h1.position - h1.normal*.01, refract);
		if (depth > 0) {
			float before = shade.x;
			shade += trace(refractRay, depth - 1) * h1.material->getRefractance();
			if (before != shade.x) {
				printf("b4 : %f\n\n", before);
				printf("af : %f\n\n", shade.x);
			}
		}

		//}


		return shade;
		//return h1.material->getColor(h1.position, h1.normal, -ray.dir);// +h1.material->getReflectionDir(-ray.dir, h1.normal);
	}

	//float3 cast(const Ray& ray) {
	//	Hit hit = firstIntersect(ray);
	//	// hit provides position r, normal n, material

	//	// if nothing hit, return background color
	//	if (hit.t < 0) return float3(1, 1, 1);


	//	return hit.material->getColor(
	//		hit.position,
	//		hit.normal,
	//		-ray.dir);
	//}

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
// global application data

// screen resolution
const int screenWidth = 600;
const int screenHeight = 600;
// image to be computed by ray tracing
float3 image[screenWidth*screenHeight];

Scene scene;

bool computeImage()
{
	static unsigned int iPart = 0;

	if(iPart >= 64)
		return false;
    for(int j = iPart; j < screenHeight; j+=64)
	{
        for(int i = 0; i < screenWidth; i++)
		{
			float3 pixelColor = float3(0, 0, 0);
			float2 ndcPixelCentre( (2.0 * i - screenWidth) / screenWidth, (2.0 * j - screenHeight) / screenHeight );

			Camera& camera = scene.getCamera();
			Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcPixelCentre));
			
			image[j*screenWidth + i] = scene.trace(ray, 5);
		}
	}
	iPart++;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL starts here. In the ray tracing example, OpenGL just outputs the image computed to the array.

// display callback invoked when window needs to be redrawn
void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear screen

	if(computeImage())
		glutPostRedisplay();
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);
 
    glutSwapBuffers(); // drawing finished
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);						// initialize GLUT
    glutInitWindowSize(screenWidth, screenHeight);				// startup window size 
    glutInitWindowPosition(100, 100);           // where to put window on screen
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    // 8 bit R,G,B,A + double buffer + depth buffer
 
    glutCreateWindow("Ray caster");				// application window is created and displayed
 
    glViewport(0, 0, screenWidth, screenHeight);

    glutDisplayFunc(onDisplay);					// register callback
 
    glutMainLoop();								// launch event handling loop
    
    return 0;
}

