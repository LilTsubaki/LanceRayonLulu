
// This code is highly based on smallpt
// http://www.kevinbeason.com/smallpt/
#include <cmath>
#include <algorithm>
#include <cassert>
#include <random>
#include <memory>
#include <fstream>
#include <iostream>
#include "maillage.h"
#include "plane.h"

// GLM (vector / matrix)
#define GLM_FORCE_RADIANS

#include <glm/vec4.hpp>
#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>

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

};
struct Box
{
    glm::vec3 v0, v1;
    Box( glm::vec3 v, glm::vec3 w):v0(v),v1(w){}

};

bool inDaBox(const Box & b, QVector3D p)
{
    return(p.x() > b.v0.x && p.x() < b.v1.x && p.y() > b.v0.y && p.y() < b.v1.y && p.z() > b.v0.z && p.z() < b.v1.z);
}

struct BoiteEnglobante
{
    BoiteEnglobante(const Maillage & m, const Box & b, const QVector<int> & indices):maillage(&m),laBoite(b)
    {
        for(int cpt = 0; cpt < indices.length(); ++cpt)
        {
            int indice=indices.at(cpt);
            if(inDaBox(b,m.getGeom().at(m.getTopo().at(3*indice)))
                    ||inDaBox(b,m.getGeom().at(m.getTopo().at(3*indice+1)))
                    ||inDaBox(b,m.getGeom().at(m.getTopo().at(3*indice+2))))
            {
                //std::cout << "append" << std::endl;
                indexTriangles.append(indice);
            }

        }
        std::cout<<indexTriangles.size()<<std::endl;

        if(indexTriangles.size()>20)
        {
            //std::cout << "subdivision" << std::endl;
            glm::vec3 tem=((b.v0+b.v1)/2.0f);
            Box b1(b.v0,tem);
            Box b2(glm::vec3(b1.v1.x, b1.v0.y, b1.v0.z),glm::vec3(b.v1.x, b1.v1.y, b1.v1.z));
            Box b3(glm::vec3(b1.v0.x, b1.v1.y, b1.v0.z),glm::vec3(b1.v1.x, b.v1.y, b1.v1.z));
            Box b4(glm::vec3(b1.v1.x, b1.v1.y, b.v0.z),glm::vec3(b.v1.x, b.v1.y, b1.v1.z));

            Box b5(glm::vec3(b.v0.x, b.v0.y, b1.v1.z),glm::vec3(b1.v1.x, b1.v1.y, b.v1.z));
            Box b6(glm::vec3(b1.v1.x, b1.v0.y, b1.v1.z),glm::vec3(b.v1.x, b1.v1.y, b.v1.z));
            Box b7(glm::vec3(b1.v0.x, b1.v1.y, b1.v1.z),glm::vec3(b1.v1.x, b.v1.y, b.v1.z));
            Box b8(glm::vec3(b1.v1.x, b1.v1.y, b1.v1.z),glm::vec3(b.v1.x, b.v1.y, b.v1.z));
            f1=new BoiteEnglobante(m,b1,indexTriangles);
            f2=new BoiteEnglobante(m,b2,indexTriangles);
            f3=new BoiteEnglobante(m,b3,indexTriangles);
            f4=new BoiteEnglobante(m,b4,indexTriangles);
            f5=new BoiteEnglobante(m,b5,indexTriangles);
            f6=new BoiteEnglobante(m,b6,indexTriangles);
            f7=new BoiteEnglobante(m,b7,indexTriangles);
            f8=new BoiteEnglobante(m,b8,indexTriangles);

        }

    }

    Box laBoite;
    BoiteEnglobante* f1;
    BoiteEnglobante* f2;
    BoiteEnglobante* f3;
    BoiteEnglobante* f4;
    BoiteEnglobante* f5;
    BoiteEnglobante* f6;
    BoiteEnglobante* f7;
    BoiteEnglobante* f8;
    const Maillage * maillage;
    QVector<int> indexTriangles;
    ~BoiteEnglobante()
    {
        delete(f1);
        delete(f2);
        delete(f3);
        delete(f4);
        delete(f5);
        delete(f6);
        delete(f7);
        delete(f8);
    }
};

struct MaillageBox
{
    const Maillage m;
    const BoiteEnglobante b;
};

glm::vec3 normale(const glm::vec3& point, const Sphere& sphere){
    return glm::normalize(point-sphere.center);
}

glm::vec3 normale(const glm::vec3& point, const Triangle &triangle){
    return glm::normalize(glm::cross(triangle.v1-triangle.v0,triangle.v2-triangle.v0));
}

glm::vec3 normale(const glm::vec3& point, const Maillage &maillage)
{
    return glm::vec3(1,0,0);
}

glm::vec3 normale(const glm::vec3& point, const BoiteEnglobante &be)
{
    return glm::vec3(1,0,0);
}

glm::vec3 normale(const glm::vec3& point, const MaillageBox &mb)
{
    return normale(point,mb.m);
}

glm::vec3 normale(const glm::vec3& point, const Box &box){
    if(std::abs(point.x-box.v0.x)<0.01)
        return glm::vec3(-1,0,0);
    if(std::abs(point.x-box.v1.x)<0.01)
        return glm::vec3(1,0,0);
    if(std::abs(point.y-box.v0.y)<0.01)
        return glm::vec3(0,-1,0);
    if(std::abs(point.y-box.v1.y)<0.01)
        return glm::vec3(0,1,0);
    if(std::abs(point.z-box.v0.z)<0.01)
        return glm::vec3(0,0,-1);
    if(std::abs(point.z-box.v1.z)<0.01)
        return glm::vec3(0,0,1);
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

float intersect (const Ray & ray, const Maillage &maillage)
{
    float t=noIntersect;
    for(int i=0;i<maillage.getTopo().size();i+=3){
        QVector3D p1=maillage.getGeom().at(maillage.getTopo().at(i));
        QVector3D p2=maillage.getGeom().at(maillage.getTopo().at(i+1));
        QVector3D p3=maillage.getGeom().at(maillage.getTopo().at(i+2));
        Triangle tri{{p1.x(),p1.y(),p1.z()},{p2.x(),p2.y(),p2.z()},{p3.x(),p3.y(),p3.z()}};
        float temp=intersect(ray,tri);
        if(temp<t&&temp>0)
            t=temp;
    }
     return t;
}

float intersect (const Ray & r, const Box &box)
{
    float t;
    glm::vec3 dirfrac;
    // r.dir is unit direction vector of ray
    dirfrac.x = 1.0f / r.direction.x;
    dirfrac.y = 1.0f / r.direction.y;
    dirfrac.z = 1.0f / r.direction.z;
    // lb is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
    // r.org is origin of ray
    float t1 = (box.v0.x - r.origin.x)*dirfrac.x;
    float t2 = (box.v1.x - r.origin.x)*dirfrac.x;
    float t3 = (box.v0.y - r.origin.y)*dirfrac.y;
    float t4 = (box.v1.y - r.origin.y)*dirfrac.y;
    float t5 = (box.v0.z - r.origin.z)*dirfrac.z;
    float t6 = (box.v1.z - r.origin.z)*dirfrac.z;

    float tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
    float tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

    // if tmax < 0, ray (line) is intersecting AABB, but whole AABB is behing us
    if (tmax < 0)
    {
        t = tmax;
        return noIntersect;
    }

    // if tmin > tmax, ray doesn't intersect AABB
    if (tmin > tmax)
    {
        t = tmax;
        return noIntersect;
    }

    t = tmin;
    return t;

}

float intersect (const Ray & ray, const BoiteEnglobante &b)
{
    float f=noIntersect;
    if(intersect(ray, b.laBoite) != noIntersect)
    {
        if(b.f1!=nullptr)
        {

            float minf;
            if((f=intersect(ray, *(b.f1)))!=noIntersect){
               if(f<minf)
                   minf=f;
            }
            if((f=intersect(ray, *(b.f2)))!=noIntersect){
                if(f<minf)
                    minf=f;
            }
            if((f=intersect(ray, *(b.f3)))!=noIntersect){
                if(f<minf)
                    minf=f;
            }
            if((f=intersect(ray, *(b.f4)))!=noIntersect){
                if(f<minf)
                    minf=f;
            }
            if((f=intersect(ray, *(b.f5)))!=noIntersect){
                if(f<minf)
                    minf=f;
            }
            if((f=intersect(ray, *(b.f6)))!=noIntersect){
                if(f<minf)
                    minf=f;
            }
            if((f=intersect(ray, *(b.f7)))!=noIntersect){
                if(f<minf)
                    minf=f;
            }
            if((f=intersect(ray, *(b.f8)))!=noIntersect){
                if(f<minf)
                    minf=f;
            }
            return minf;
        }
        else
        {
            for(int cpt = 0; cpt < b.indexTriangles.length(); cpt++)
            {
                QVector3D p1=b.maillage->getGeom().at(b.maillage->getTopo().at(b.indexTriangles.at(cpt)*3));
                QVector3D p2=b.maillage->getGeom().at(b.maillage->getTopo().at(b.indexTriangles.at(cpt)*3+1));
                QVector3D p3=b.maillage->getGeom().at(b.maillage->getTopo().at(b.indexTriangles.at(cpt)*3+2));
                Triangle tri{{p1.x(),p1.y(),p1.z()},{p2.x(),p2.y(),p2.z()},{p3.x(),p3.y(),p3.z()}};
                float temp=intersect(ray,tri);
                if(temp<f&&temp>0)
                    f=temp;
            }

            return f;
        }
    }
    else
    {
        return noIntersect;
    }
}

float intersect (const Ray & ray, const MaillageBox &mb)
{
    return intersect(ray,mb.b);
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
glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Diffuse& mat,int it);
glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Glass& mat,int it);
glm::vec3 directLight(const Ray& r, const glm::vec3& norm, float f, const Mirror& mat, int it);
glm::vec3 reflect(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Diffuse& mat, float &fres);
glm::vec3 reflect(const Ray& r,const glm::vec3& norm,glm::vec3 &dir, const Glass& mat, float &fres);
glm::vec3 reflect(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Mirror& mat, float &fres);
glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Diffuse& mat, float &fres);
glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Glass& mat, float &fres);
glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Mirror& mat, float &fres);


struct Object
{
    virtual float intersect(const Ray &r) const = 0;
    virtual glm::vec3 albedo() const = 0;
    virtual glm::vec3 normale(const glm::vec3& point) const =0;
    virtual glm::vec3 directLight(const Ray& ray, const glm::vec3& normale, float f,int it) const =0;
    virtual glm::vec3 reflect(const Ray& ray, const glm::vec3& normale,glm::vec3 &dir, float &fres) const =0;
    virtual glm::vec3 refract(const Ray& ray, const glm::vec3& normale,glm::vec3 &dir, float &fres) const =0;
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

    glm::vec3 directLight(const Ray& ray, const glm::vec3& normale, float f,int it) const
    {
        return ::directLight(ray, normale, f, material,it);
    }
    glm::vec3 reflect(const Ray& ray, const glm::vec3& normale,glm::vec3 &dir, float &fres) const{
        return ::reflect(ray, normale,dir, material, fres);
    }

    glm::vec3 refract(const Ray& ray, const glm::vec3& normale,glm::vec3 &dir, float &fres) const{
        return ::refract(ray, normale,dir,material,fres);
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
    const Triangle leftWallA{{0, 0, 0}, {0, 100, 0},{0, 0, 150}};
    const Triangle leftWallB{{0, 100, 150}, {0, 0, 150}, {0, 100, 0}};

    // Right Wall
    const Triangle rightWallA{{100, 0, 0}, {100, 0, 150}, {100, 100, 0}};
    const Triangle rightWallB{{100, 100, 150}, {100, 100, 0}, {100, 0, 150}};

    // Back wall
    const Triangle backWallA{{0, 0, 0}, {100, 0, 0}, {100, 100, 0}};
    const Triangle backWallB{{0, 0, 0}, {100, 100, 0}, {0, 100, 0}};

    // Bottom Floor
    const Triangle bottomWallA{{0, 0, 0}, {100, 0, 150}, {100, 0, 0}};
    const Triangle bottomWallB{{0, 0, 0}, {0, 0, 150}, {100, 0, 150}};

    // Top Ceiling
    const Triangle topWallA{{0, 100, 0}, {100, 100, 0}, {0, 100, 150}};
    const Triangle topWallB{{100, 100, 150}, {0, 100, 150}, {100, 100, 0}};

    const Sphere leftSphere{16.5, glm::vec3 {27, 16.5, 47}};
    const Sphere rightSphere{16.5, glm::vec3 {73, 16.5, 78}};
    //const Sphere sphere{25.0, glm::vec3 {27, 45, 47}};

    //const glm::vec3 light{27, 45, 47};
    const glm::vec3 light{50, 70, 81.6};
    const glm::vec3 lightColor{1,1,1};

    // Materials
    const Diffuse white{{.75, .75, .75}};
    const Diffuse red{{.75, .25, .25}};
    const Diffuse blue{{.25, .25, .75}};

    const Glass glass{{1, 1, 1}};
    const Mirror mirror{{1, 1, 1}};

    //const Box box{glm::vec3 {50, 16.5, 78}, glm::vec3 {55, 21.5, 83}};



    // Objects
    // Note: this is a rather convoluted way of initialising a vector of unique_ptr ;)
    const std::vector<std::unique_ptr<Object>> objects = [] (){
        Maillage m;
        QVector3D center(50,0,50);
        glm::vec3 min;
        glm::vec3 max;
        m.geometry(center,"C:/Users/etu/Downloads/BeautifulGirl.obj", min, max);
        std::cout << "chargement obj done " << std::endl;
        m.translate(center, min, max);
        //m.Ecriture("BeautifulGirl.obj");
        Box box(min, max);
        std::cout << "debut creation boite englobante " << std::endl;
        BoiteEnglobante test(m, box, m.getTopo());
        std::cout << test.indexTriangles.length() << std::endl;
        std::cout << test.f1->indexTriangles.length() << std::endl;
        std::cout << test.f2->indexTriangles.length() << std::endl;
        std::cout << test.f3->indexTriangles.length() << std::endl;
        std::cout << test.f4->indexTriangles.length() << std::endl;
        std::cout << test.f5->indexTriangles.length() << std::endl;
        std::cout << test.f6->indexTriangles.length() << std::endl;
        std::cout << test.f7->indexTriangles.length() << std::endl;
        std::cout << test.f8->indexTriangles.length() << std::endl;
        std::cout << "creation boite englobante done " << std::endl;

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

        //ret.push_back(makeObject(m, red));

        //ret.push_back(makeObject(m, white));
        //ret.push_back(std::unique_ptr<Object>(new ObjectTpl<Maillage,Diffuse>(m,white)));

        /*for(int i=0;i<m.getTopo().size();i+=3){
            QVector3D p1=m.getGeom().at(m.getTopo().at(i));
            QVector3D p2=m.getGeom().at(m.getTopo().at(i+1));
            QVector3D p3=m.getGeom().at(m.getTopo().at(i+2));
            Triangle tri{{p1.x(),p1.y(),p1.z()},{p2.x(),p2.y(),p2.z()},{p3.x(),p3.y(),p3.z()}};
            ret.push_back(makeObject(tri,white));
        }*/
        //ret.push_back(makeObject(leftSphere, mirror));
        //ret.push_back(makeObject(rightSphere, glass));
        //ret.push_back(makeObject(sphere, white));

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

glm::vec3 radiance (const Ray & r,int i=8)
{
    float fres;
    float f;
    glm::vec3 norm;
    glm::vec3 dirReflect;
    glm::vec3 dirRefract;

    Object* obj=intersect(r,f,norm);
    if(obj == nullptr)
        return glm::vec3(0,0,0);
    glm::vec3 color=obj->directLight(r,norm,f,i);

    glm::vec3 colorReflect(0,0,0);
    glm::vec3 colorRefract(0,0,0);
    if(i>0){
        colorReflect=obj->reflect(r,norm,dirReflect, fres);

        if(random_u()<fres){

            glm::vec3 ori=r.origin+f*r.direction+0.1f*dirReflect;
            Ray ray{ori,dirReflect};

            colorReflect=colorReflect*radiance(ray,i-1);
            colorRefract=glm::vec3(0,0,0);
        }
        else{
            colorRefract=obj->refract(r,norm,dirRefract, fres);
            glm::vec3 orig=r.origin+f*r.direction+0.1f*dirRefract;
            Ray ray2{orig,dirRefract};

            colorRefract=colorRefract*radiance(ray2,i-1);
            colorReflect=glm::vec3(0,0,0);
        }
    }

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


    int cmpt=0;

    #pragma omp parallel for //schedule(dynamic,5)
    for (int y = 0; y < h; y++)
    {
        std::cerr << "\rRendering: " << 100 * cmpt++ / (h - 1) << "%";

        for (unsigned short x = 0; x < w; x++)
        {
            int al=8;
            glm::vec3 r(0,0,0);
            for(int i=0;i<al;i++)
            {
                float random1 = random_u();
                float random2 = random_u();
            glm::vec4 p0 = screenToRay * glm::vec4{float(x+random1-0.5), float(h - y+random2-0.5), 0.f, 1.f};
            glm::vec4 p1 = screenToRay * glm::vec4{float(x+random1-0.5), float(h - y+random2-0.5), 1.f, 1.f};

            glm::vec3 pp0 = glm::vec3(p0 / p0.w);
            glm::vec3 pp1 = glm::vec3(p1 / p1.w);

            glm::vec3 d = glm::normalize(pp1 - pp0);
            r += radiance (Ray{pp0, d});
            }
             r/=al;
            colors[y * w + x] += glm::clamp(r, glm::vec3(0.f, 0.f, 0.f), glm::vec3(1.f, 1.f, 1.f));
        }
    }

    {
        std::fstream f("image.ppm", std::fstream::out);
        f << "P3\n" << w << " " << h << std::endl << "255" << std::endl;

        for (auto c : colors)
            f << toInt(c.x) << " " << toInt(c.y) << " " << toInt(c.z) << " ";
    }
}

glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Diffuse& mat,int it){
    glm::vec3 norm2;
    glm::vec3 color=mat.color;
    glm::vec3 inter=r.origin+f*r.direction;
    float pdf;

    glm::vec3 dirInterToLight=scene::light-inter;
    glm::vec3 dirInterToPointLight=glm::normalize(dirInterToLight);
    dirInterToPointLight=scene::light+sample_sphere(10.0f,random_u(),random_u(),pdf,-dirInterToPointLight)-inter;
    glm::vec3 dirInterToPointLightNormal=glm::normalize(dirInterToPointLight);
    Ray rayInterToPointLight{inter+(0.1f*dirInterToPointLightNormal),dirInterToPointLightNormal};
    float distInterToObj;
    intersect(rayInterToPointLight,distInterToObj,norm2);
    if(distInterToObj*distInterToObj<glm::dot(dirInterToPointLight,dirInterToPointLight)){
                color=glm::vec3(0,0,0);
    }



    float fact =std::abs(glm::dot(dirInterToPointLightNormal,norm)/pi);


    return (fact*color/pdf/(glm::dot(dirInterToPointLight,dirInterToPointLight)))*scene::lightColor;
}

glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Glass& mat,int it){

    return glm::vec3(0,0,0);
}

glm::vec3 directLight(const Ray& r,const glm::vec3& norm, float f, const Mirror& mat,int it){
    return glm::vec3(0,0,0);
}

glm::vec3 reflect(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Diffuse& mat, float &fres){
    float pdf2;
    dir=glm::normalize(sample_sphere(1,random_u(),random_u(),pdf2, norm));
    fres=1;
    return mat.color;
}

glm::vec3 reflect(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Glass& mat, float &fres){
    dir=reflect(-r.direction,norm);
    fres=fresnelR(-r.direction,norm,1.5f);
    return mat.color;//*/
    //return glm::vec3(0,0,0);
}

glm::vec3 reflect(const Ray& r, const glm::vec3& norm, glm::vec3& dir, const Mirror& mat, float &fres){
       fres = 1;
    dir=reflect(-r.direction,norm);
    return mat.color;
}

glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Diffuse& mat, float &fres){
    fres = 1;
    return glm::vec3(0,0,0);
}

glm::vec3 refract(const Ray& r,const glm::vec3& norm,glm::vec3 &dir, const Glass& mat, float &fres){
    if(refract(-r.direction,norm,1.5f,dir)){
        fres=fresnelR(-r.direction,norm,1.5f);
        return mat.color;
    }
    else{
        fres = 1;
        return glm::vec3(0,0,0);
    }//*/

    //return glm::vec3(0,0,0);
}

glm::vec3 refract(const Ray& r, const glm::vec3& norm, glm::vec3 &dir, const Mirror& mat, float &fres){
    fres = 1;
    return glm::vec3(0,0,0);
}

