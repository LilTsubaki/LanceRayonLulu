// This code is highly based on smallpt
// http://www.kevinbeason.com/smallpt/
#include <cmath>
#include <algorithm>
#include <cassert>
#include <random>
#include <memory>
#include <fstream>
#include <iostream>

// GLM (vector / matrix)
#define GLM_FORCE_RADIANS

#include <C:/Users/etu/Downloads/glm-0.9.7.1/glm/glm/vec4.hpp>
#include <C:/Users/etu/Downloads/glm-0.9.7.1/glm/glm/vec3.hpp>
#include <C:/Users/etu/Downloads/glm-0.9.7.1/glm/glm/mat4x4.hpp>
#include <C:/Users/etu/Downloads/glm-0.9.7.1/glm/glm/gtc/matrix_transform.hpp>

const float pi = 3.1415927f;
const float noIntersect = std::numeric_limits<float>::infinity();

bool isIntersect(float t)
{
	return t < noIntersect;
}

struct Ray
{
	const glm::vec3 origin, direction;
};

struct Sphere
{
	const float radius;
	const glm::vec3 center;
};

struct Triangle
{
	const glm::vec3 v0, v1, v2;
        //Triangle(glm::vec3 va,glm::vec3 vb, glm::vec3 vc ) : v0(va), v1(vb), v2(vc){}

};

glm::vec3 normale(const glm::vec3& point, const Sphere& sphere){
    return glm::normalize(point-sphere.center);
}

glm::vec3 normale(const glm::vec3& point, const Triangle &triangle){
    return glm::normalize(glm::cross(triangle.v1-triangle.v0,triangle.v2-triangle.v0));
}

	// WARRING: works only if r.d is normalized
float intersect (const Ray & ray, const Sphere &sphere)
{				// returns distance, 0 if nohit
	glm::vec3 op = sphere.center - ray.origin;		// Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
	float t, b = glm:: dot(ray.direction, op), det =
		b * b - glm::dot(op, op) + sphere.radius * sphere.radius;
	if (det < 0)
		return noIntersect;
	else
		det = std::sqrt (det);
	return (t = b - det) >= 0 ? t : ((t = b + det) >= 0 ? t : noIntersect);
}

float intersect(const Ray & ray, const Triangle &triangle)
{
	auto e1 = triangle.v1 - triangle.v0;
	auto e2 = triangle.v2 - triangle.v0;

	auto h = glm::cross(ray.direction, e2);
	auto a = glm::dot(e1, h);

	auto f = 1.f / a;
	auto s = ray.origin - triangle.v0;

	auto u = f * glm::dot(s, h);
	auto q = glm::cross(s, e1);
	auto v = f * glm::dot(ray.direction, q);
	auto t = f * glm::dot(e2, q);

	if(std::abs(a) < 0.00001)
		return noIntersect;
	if(u < 0 || u > 1)
		return noIntersect;
	if(v < 0 || (u+v) > 1)
		return noIntersect;
	if(t < 0)
		return noIntersect;

	return t;
}

struct Diffuse
{
	const glm::vec3 color;
};

struct Glass
{
	const glm::vec3 color;
};

struct Mirror
{
	const glm::vec3 color;
};

template<typename T>
glm::vec3 albedo(const T &t)
{
	return t.color;
}
glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Diffuse& mat);
glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Glass& mat);
glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Mirror& mat);
glm::vec3 reflect(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Diffuse& mat);
glm::vec3 reflect(const Ray& r,const glm::vec3& norm,glm::vec3 &dir, const Glass& mat);
glm::vec3 reflect(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Mirror& mat);
glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Diffuse& mat);
glm::vec3 refract(const Ray& r,const glm::vec3& norm,glm::vec3 &dir, const Glass& mat);
glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Mirror& mat);


struct Object
{
	virtual float intersect(const Ray &r) const = 0;
	virtual glm::vec3 albedo() const = 0;
    virtual glm::vec3 normale(const glm::vec3& point) const =0;
    virtual glm::vec3 directLight(const Ray& ray, const glm::vec3& normale, float f) const =0;
    virtual glm::vec3 reflect(const Ray& ray, const glm::vec3& normale,glm::vec3 &dir) const =0;
    virtual glm::vec3 refract(const Ray& ray, const glm::vec3& normale,glm::vec3 &dir) const =0;
};

template<typename P, typename M>
struct ObjectTpl final : Object
{
	ObjectTpl(const P &_p, const M &_m)
		:primitive(_p), material(_m)
	{}

	float intersect(const Ray &ray) const
	{
		return ::intersect(ray, primitive);
	}

    glm::vec3 normale(const glm::vec3& point) const
    {
        return ::normale(point, primitive);
    }

	glm::vec3 albedo() const
	{
		return ::albedo(material);
	}

    glm::vec3 directLight(const Ray& ray, const glm::vec3& normale, float f) const
    {
        return ::directLight(ray, normale, f, material);
    }
    glm::vec3 reflect(const Ray& ray, const glm::vec3& normale,glm::vec3 &dir) const{
        return ::reflect(ray, normale,dir, material);
    }

    glm::vec3 refract(const Ray& ray, const glm::vec3& normale,glm::vec3 &dir) const{
        return ::refract(ray, normale,dir,material);
    }

	const P &primitive;
	const M &material;
};


template<typename P, typename M>
std::unique_ptr<Object> makeObject(const P&p, const M&m)
{
	return std::unique_ptr<Object>(new ObjectTpl<P, M>{p, m});
}

// Scene
namespace scene
{
	// Primitives

	// Left Wall
    const Triangle leftWallA{{0, 0, 0}, {0, 100, 0}, {0, 0, 150}};
	const Triangle leftWallB{{0, 100, 150}, {0, 100, 0}, {0, 0, 150}};

	// Right Wall
	const Triangle rightWallA{{100, 0, 0}, {100, 100, 0}, {100, 0, 150}};
	const Triangle rightWallB{{100, 100, 150}, {100, 100, 0}, {100, 0, 150}};

	// Back wall
	const Triangle backWallA{{0, 0, 0}, {100, 0, 0}, {100, 100, 0}};
	const Triangle backWallB{{0, 0, 0}, {0, 100, 0}, {100, 100, 0}};

	// Bottom Floor
	const Triangle bottomWallA{{0, 0, 0}, {100, 0, 0}, {100, 0, 150}};
	const Triangle bottomWallB{{0, 0, 0}, {0, 0, 150}, {100, 0, 150}};

	// Top Ceiling
	const Triangle topWallA{{0, 100, 0}, {100, 100, 0}, {0, 100, 150}};
	const Triangle topWallB{{100, 100, 150}, {100, 100, 0}, {0, 100, 150}};

	const Sphere leftSphere{16.5, glm::vec3 {27, 16.5, 47}};
	const Sphere rightSphere{16.5, glm::vec3 {73, 16.5, 78}};

	const glm::vec3 light{50, 70, 81.6};

	// Materials
	const Diffuse white{{.75, .75, .75}};
	const Diffuse red{{.75, .25, .25}};
	const Diffuse blue{{.25, .25, .75}};

    const Glass glass{{1, 1, 1}};
    const Mirror mirror{{1, 1, 1}};

	// Objects
	// Note: this is a rather convoluted way of initialising a vector of unique_ptr ;)
	const std::vector<std::unique_ptr<Object>> objects = [] (){
		std::vector<std::unique_ptr<Object>> ret;
		ret.push_back(makeObject(backWallA, white));
		ret.push_back(makeObject(backWallB, white));
		ret.push_back(makeObject(topWallA, white));
		ret.push_back(makeObject(topWallB, white));
		ret.push_back(makeObject(bottomWallA, white));
		ret.push_back(makeObject(bottomWallB, white));
		ret.push_back(makeObject(rightWallA, blue));
		ret.push_back(makeObject(rightWallB, blue));
		ret.push_back(makeObject(leftWallA, red));
		ret.push_back(makeObject(leftWallB, red));

		ret.push_back(makeObject(leftSphere, mirror));
        ret.push_back(makeObject(rightSphere, glass));

		return ret;
	}();
}

thread_local std::default_random_engine generator;
thread_local std::uniform_real_distribution<float> distribution(0.0,1.0);

float random_u()
{
	return distribution(generator);
}

glm::vec3 sample_cos(const float u, const float v, const glm::vec3 n)
{
	// Ugly: create an ornthogonal base
	glm::vec3 basex, basey, basez;

	basez = n;
	basey = glm::vec3(n.y, n.z, n.x);

	basex = glm::cross(basez, basey);
	basex = glm::normalize(basex);

	basey = glm::cross(basez, basex);

	// cosinus sampling. Pdf = cosinus
	return  basex * (std::cos(2.f * pi * u) * std::sqrt(1.f - v)) +
		basey * (std::sin(2.f * pi * u) * std::sqrt(1.f - v)) +
		basez * std::sqrt(v);
}

int toInt (const float x)
{
	return int (std::pow (glm::clamp (x, 0.f, 1.f), 1.f / 2.2f) * 255 + .5);
}

// WARNING: ASSUME NORMALIZED RAY
// Compute the intersection ray / scene.
// Returns true if intersection
// t is defined as the abscisce along the ray (i.e
//             p = r.o + t * r.d
// id is the id of the intersected object
Object* intersect (const Ray & r, float &t, glm::vec3& norm)
{
	t = noIntersect;
	Object *ret = nullptr;

	for(auto &object : scene::objects)
	{
		float d = object->intersect(r);
		if (isIntersect(d) && d < t)
		{
			t = d;
			ret = object.get();
		}
	}
    glm::vec3 point(r.origin+t*r.direction);
    if(ret!=nullptr)
        norm=ret->normale(point);

	return ret;
}

// Reflect the ray i along the normal.
// i should be oriented as "leaving the surface"
glm::vec3 reflect(const glm::vec3 i, const glm::vec3 n)
{
	return n * (glm::dot(n, i)) * 2.f - i;
}

float sin2cos (const float x)
{
	return std::sqrt(std::max(0.0f, 1.0f-x*x));
}

// Fresnel coeficient of transmission.
// Normal point outside the surface
// ior is n0 / n1 where n0 is inside and n1 is outside
float fresnelR(const glm::vec3 i, const glm::vec3 n, const float ior)
{
	if(glm::dot(n, i) < 0)
		return fresnelR(i, n * -1.f, 1.f / ior);

	float R0 = (ior - 1.f) / (ior + 1.f);
	R0 *= R0;

	return R0 + (1.f - R0) * std::pow(1.f - glm::dot(i, n), 5.f);
}

// compute refraction vector.
// return true if refraction is possible.
// i and n are normalized
// output wo, the refracted vector (normalized)
// n point oitside the surface.
// ior is n00 / n1 where n0 is inside and n1 is outside
//
// i point outside of the surface
bool refract(glm::vec3 i, glm::vec3 n, float ior, glm::vec3 &wo)
{
	i = i * -1.f;

	if(glm::dot(n, i) > 0)
	{
		n = n * -1.f;
	}
	else
	{
		ior = 1.f / ior;
	}

    float k = 1.f - ior * ior * (1.f - glm::dot(n, i) * glm::dot(n, i));
    if (k < 0.)
		return false;

	wo = i * ior - n * (ior * glm::dot(n, i) + std::sqrt(k));

	return true;
}

glm::vec3 sample_sphere(const float r, const float u, const float v, float &pdf, const glm::vec3& normal)
{
	pdf = 1.f / (pi * r * r);
	glm::vec3 sample_p = sample_cos(u, v, normal);

	float cos = glm::dot(sample_p, normal);

	pdf *= cos;
	return sample_p * r;
}

glm::vec3 radiance (const Ray & r,int i=0)
{
    float f;
    glm::vec3 norm;
    glm::vec3 dir;

    Object* obj=intersect(r,f,norm);
    if(obj == nullptr)
        return glm::vec3(0,0,0);
    glm::vec3 color=obj->directLight(r,norm,f);
    glm::vec3 colorReflect=obj->reflect(r,norm,dir);
    glm::vec3 ori=r.origin+f*r.direction+0.1f*dir;
    Ray ray{ori,dir};
    if(i<8)//||colorReflect.length()<0.01)
        colorReflect=colorReflect*radiance(ray,i+1);
    else
        colorReflect=glm::vec3(0,0,0);


    glm::vec3 colorRefract=obj->refract(r,norm,dir);
    glm::vec3 orig=r.origin+f*r.direction+0.1f*dir;
    Ray ray2{orig,dir};
    if(i<8)//||colorRefract.length()<0.01)
        colorRefract=colorRefract*radiance(ray2,i+1);
    else
        colorRefract=glm::vec3(0,0,0);

    //glm::vec3 color = obj->albedo();
    /*glm::vec3 inter=r.origin+f*r.direction;
    glm::vec3 dir=scene::light-inter;
    glm::vec3 dirn=glm::normalize(dir);
    Ray s{inter+(0.1f*dirn),dirn};
    intersect(s,f,norm2);
    if(f*f<glm::dot(dir,dir)){
        return color*0.1f;
        //return glm::vec3(0,0,0);
    }
    float fact =glm::dot(dirn,norm)/pi;


    if(fact<0) fact=-fact;
    fact*=0.9;
    fact+=0.1;

    return fact*color;*/
    return color+colorReflect+colorRefract;
}

int main (int, char **)
{
	int w = 768, h = 768;
	std::vector<glm::vec3> colors(w * h, glm::vec3{0.f, 0.f, 0.f});

	Ray cam {{50, 52, 295.6}, glm::normalize(glm::vec3{0, -0.042612, -1})};	// cam pos, dir
	float near = 1.f;
	float far = 10000.f;

	glm::mat4 camera =
		glm::scale(glm::mat4(1.f), glm::vec3(float(w), float(h), 1.f))
		* glm::translate(glm::mat4(1.f), glm::vec3(0.5, 0.5, 0.f))
		* glm::perspective(float(54.5f * pi / 180.f), float(w) / float(h), near, far)
		* glm::lookAt(cam.origin, cam.origin + cam.direction, glm::vec3(0, 1, 0))
		;

	glm::mat4 screenToRay = glm::inverse(camera);

	for (int y = 0; y < h; y++)
    {
		std::cerr << "\rRendering: " << 100 * y / (h - 1) << "%";

		for (unsigned short x = 0; x < w; x++)
		{
			glm::vec4 p0 = screenToRay * glm::vec4{float(x), float(h - y), 0.f, 1.f};
			glm::vec4 p1 = screenToRay * glm::vec4{float(x), float(h - y), 1.f, 1.f};

			glm::vec3 pp0 = glm::vec3(p0 / p0.w);
			glm::vec3 pp1 = glm::vec3(p1 / p1.w);

			glm::vec3 d = glm::normalize(pp1 - pp0);

			glm::vec3 r = radiance (Ray{pp0, d});
			colors[y * w + x] += glm::clamp(r, glm::vec3(0.f, 0.f, 0.f), glm::vec3(1.f, 1.f, 1.f)) * 0.25f;
		}
    }

	{
		std::fstream f("image.ppm", std::fstream::out);
		f << "P3\n" << w << " " << h << std::endl << "255" << std::endl;

		for (auto c : colors)
			f << toInt(c.x) << " " << toInt(c.y) << " " << toInt(c.z) << " ";
	}
}

glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Diffuse& mat){
    glm::vec3 norm2;
    glm::vec3 color=mat.color;
    glm::vec3 inter=r.origin+f*r.direction;
    //random
    glm::vec3 dir=scene::light-inter;
    glm::vec3 dirn=glm::normalize(dir);
    Ray s{inter+(0.1f*dirn),dirn};
    intersect(s,f,norm2);
    if(f*f<glm::dot(dir,dir)){
        return glm::vec3(0,0,0);
    }
    float fact =std::abs(glm::dot(dirn,norm)/pi);


    return fact*color;
}

glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Glass& mat){
    /*glm::vec3 norm2;
    glm::vec3 color=mat.color;
    glm::vec3 inter=r.origin+f*r.direction;
    glm::vec3 dir=scene::light-inter;
    glm::vec3 dirn=glm::normalize(dir);
    Ray s{inter+(0.1f*dirn),dirn};
    intersect(s,f,norm2);
    if(f*f<glm::dot(dir,dir)){
        return glm::vec3(0,0,0);
    }
    float fact =std::abs(glm::dot(dirn,norm)/pi);


    return fact*color;//*/
    return glm::vec3(0,0,0);
}

glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Mirror& mat){
    return glm::vec3(0,0,0);
}

glm::vec3 reflect(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Diffuse& mat){
    return glm::vec3(0,0,0);
}

glm::vec3 reflect(const Ray& r,const glm::vec3& norm,glm::vec3 &dir, const Glass& mat){
    dir=reflect(-r.direction,norm);
    float fres=fresnelR(-r.direction,norm,1.5f);
    return fres*mat.color;//*/
    //return glm::vec3(0,0,0);
}

glm::vec3 reflect(const Ray& r,const glm::vec3& norm,glm::vec3& dir, const Mirror& mat){
    dir=reflect(-r.direction,norm);
    return mat.color;
}

glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Diffuse& mat){
    return glm::vec3(0,0,0);
}

glm::vec3 refract(const Ray& r,const glm::vec3& norm,glm::vec3 &dir, const Glass& mat){
    if(refract(-r.direction,norm,1.5f,dir)){
        float fres=fresnelR(-r.direction,norm,1.5f);
        return (1-fres)*mat.color;
    }
    else{
        return glm::vec3(0,0,0);
    }//*/

    //return glm::vec3(0,0,0);
}

glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Mirror& mat){
    return glm::vec3(0,0,0);
}
