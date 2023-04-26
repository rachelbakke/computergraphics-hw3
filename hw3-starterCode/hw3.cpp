/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Rachel Bakke
 * *************************
 * ./hw3 ./spheres.scene output.jpg
 *  ./hw3 ./test1.scene output.jpg
 * ./hw3 ./test2.scene output.jpg
 * ./hw3 ./table.scene output.jpg
 * ./hw3 ./SIGGRAPH.scene output.jpg
 * ./hw3 ./toy.scene output.jpg
 */

#ifdef WIN32
#include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
#include <GL/gl.h>
#include <GL/glut.h>
#elif defined(__APPLE__)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif
#include <glm/glm.hpp>
#include "ray.cpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <vector>
#include <string>
#include <algorithm>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char *filename = NULL;

// different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

// the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

using namespace std;

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
string intersectObjType = "sphere";
int intersectObjIndex = 0;
// glm::vec3 intersectObjPos(0, 0, 0);
bool debugMode = false;

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);
glm::vec3 colorRay(Ray curr);
float intersect(Ray curr, float lightT);

glm::vec3 maxVal(MAXFLOAT, MAXFLOAT, MAXFLOAT);
// MODIFY THIS FUNCTION
void draw_scene()
{
  float aVal = (WIDTH * 1.0f) / (HEIGHT * 1.0f);
  float topLeftX = -aVal * tan(glm::radians(fov) / 2.0f);
  float topLeftY = tan(glm::radians(fov) / 2.0f);
  float topLeftZ = -1;
  float topRightX = aVal * tan(glm::radians(fov) / 2.0f);
  float topRightY = tan(glm::radians(fov) / 2.0f);
  float topRightZ = -1;
  float botRightX = aVal * tan(glm::radians(fov) / 2.0f);
  float botRightY = -tan(glm::radians(fov) / 2.0f);
  float botRightZ = -1;
  float botLeftX = -aVal * tan(glm::radians(fov) / 2.0f);
  float botLeftY = -tan(glm::radians(fov) / 2.0f);
  float botLeftZ = -1;

  float width3D = sqrt(pow(topRightX - topLeftX, 2) + pow(topRightY - topLeftY, 2) + pow(topRightZ - topLeftZ, 2)) * 1.0f;
  float height3D = sqrt(pow(topRightX - botRightX, 2) + pow(topRightY - botRightY, 2) + pow(topRightZ - botRightZ, 2) * 1.0f);
  float wStepSize = (width3D * 1.0f) / (1.0f * WIDTH);
  float hStepSize = (height3D * 1.0f) / (HEIGHT * 1.0f);

  // get screen setup
  Ray rays[WIDTH][HEIGHT];
  for (int i = 0; i < WIDTH; i++)
  {
    for (int j = 0; j < HEIGHT; j++)
    {
      // STEP 1! generate rays
      float xVal = (i * wStepSize) + botLeftX;
      float yVal = (j * hStepSize) + botLeftY;
      glm::vec3 direction(xVal, yVal, -1.0); //
      direction = normalize(direction);
      glm::vec3 origin1((0, 0, 0));
      Ray currRay = Ray(origin1, direction);
      rays[i][j] = currRay;
    }
  }
  // PLOTTING PIXELS - get color of the pixel from the ray
  for (unsigned int x = 0; x < WIDTH; x++)
  {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for (unsigned int y = 0; y < HEIGHT; y++)
    {
      // plot_pixel(x, y, x % 256, y % 256, (x + y) % 256); //OG
      glm::vec3 currColor = colorRay(rays[x][y]);
      plot_pixel(x, y, currColor.x * 255, currColor.y * 255, currColor.z * 255); // FINAL COLOR

      // anti-aliasing
      /** glm::vec3 currColor2 = colorRay(rays[x + 1][y]);
      glm::vec3 currColor3 = colorRay(rays[x][y + 1]);
      glm::vec3 currColor4 = colorRay(rays[x + 1][y + 1]);
      float avgx = (currColor.x + currColor2.x + currColor3.x + currColor4.x) * 255 / 4;
      float avgy = (currColor.y + currColor2.y + currColor3.y + currColor4.y) * 255 / 4;
      float avgz = (currColor.z + currColor2.z + currColor3.z + currColor4.z) * 255 / 4;

      plot_pixel(x, y, avgx, avgy, avgz);
      plot_pixel(x + 1, y, avgx, avgy, avgz);
      plot_pixel(x, y + 1, avgx, avgy, avgz);
      plot_pixel(x + 1, y + 1, avgx, avgy, avgz); */

      // just intersextion
      /** glm::vec3 intersectPoint = intersect(rays[x][y]);
       if (intersectPoint != maxVal) // yes intersects
       {
         plot_pixel(x, y, 255, 192, 203); // COLOR INTERSECT
       }
       else
       {
         plot_pixel(x, y, 0.0f, 0.0f, 0.0f);
       }*/
    }
    glEnd();
    glFlush();
  }

  printf("Done!\n");
  fflush(stdout);
}

glm::vec3 colorRay(Ray curr)
{
  glm::vec3 color; // set as black
  float currT = intersect(curr, MAXFLOAT);
  glm::vec3 posObjIntersect(curr.origin.x + currT * curr.direction.x, curr.origin.y + currT * curr.direction.y, curr.origin.z + currT * curr.direction.z);
  bool checkingTri = false;
  // printf(" TEST %f , %f , %f", posIntersect.x, posIntersect.y, posIntersect.z);
  if (currT != MAXFLOAT) // object intersects
  {
    // check shadow rays for each light source and see if they are hit
    for (int l = 0; l < num_lights; l++)
    {
      // create more lights, loop through those lights,
      //  printf("\nlight: %f, %f, %f ", lights[l].position[0], lights[l].position[1], lights[l].position[2]);
      Ray currShadow = Ray();
      currShadow.setOrigin(posObjIntersect);
      glm::vec3 objToLight(lights[l].position[0] - posObjIntersect.x, lights[l].position[1] - posObjIntersect.y, lights[l].position[2] - posObjIntersect.z);
      float distanceToLight = glm::length(objToLight);
      objToLight = normalize(objToLight);

      glm::vec3 newShadowOrigin(posObjIntersect.x, posObjIntersect.y, posObjIntersect.z + 1e-5 * objToLight.z);

      currShadow.setOrigin(newShadowOrigin);
      currShadow.setDirection(objToLight);

      float shadowT = intersect(currShadow, distanceToLight);
      glm::vec3 posShadowIntersect(currShadow.origin.x + shadowT * currShadow.direction.x,
                                   currShadow.origin.y + shadowT * currShadow.direction.y,
                                   currShadow.origin.z + shadowT * currShadow.direction.z);

      glm::vec3 objToShadow(posShadowIntersect.x - posObjIntersect.x,
                            posShadowIntersect.y - posObjIntersect.y,
                            posShadowIntersect.z - posObjIntersect.z);
      float distShadow = glm::length(objToShadow);

      if (shadowT == MAXFLOAT && shadowT > 1e-3) // NOT IN SHADOW, distShadow < 1e-3
      {

        // phong model solve with normals
        // get normal of the intersection and values of shininess and everything
        if (intersectObjType == "sphere")
        {
          Sphere currSphere = spheres[intersectObjIndex];
          glm::vec3 normal(posObjIntersect.x - currSphere.position[0], posObjIntersect.y - currSphere.position[1], posObjIntersect.z - currSphere.position[2]);
          normal = glm::normalize(normal / float(currSphere.radius));
          glm::vec3 R = (2.0f * dot(objToLight, normal) * normal) - objToLight;
          R = normalize(R);
          //     clamp dot products to 0
          float kdx = float(currSphere.color_diffuse[0]);
          float kdy = float(currSphere.color_diffuse[1]);
          float kdz = float(currSphere.color_diffuse[2]);

          float ksx = float(currSphere.color_specular[0]);
          float ksy = float(currSphere.color_specular[1]);
          float ksz = float(currSphere.color_specular[2]);

          glm::vec3 v = -curr.direction;

          float RdotV = glm::dot(v, R);
          float LdotN = glm::dot(objToLight, normal);
          if (LdotN < 0)
            LdotN = 0.0;
          else if (LdotN > 1.0)
            LdotN = 1.0;
          if (RdotV < 0)
            RdotV = 0.0;
          else if (RdotV > 1.0)
            RdotV = 1.0;

          float RVdotExpo = pow(RdotV, currSphere.shininess);

          float Ix = lights[l].color[0] * (kdx * LdotN + ksx * RVdotExpo);
          float Iy = lights[l].color[1] * (kdy * LdotN + ksy * RVdotExpo);
          float Iz = lights[l].color[2] * (kdz * LdotN + ksz * RVdotExpo);
          color.x += Ix;
          color.y += Iy;
          color.z += Iz;
          // printf("\nthis is color values: %f , %f, %f ", color.x, color.y, color.z);
        }
        else
        {

          Triangle currTri = triangles[intersectObjIndex];
          //  need to know shape type and index to access it here and get its normal
          glm::vec3 e1(currTri.v[1].position[0] - currTri.v[0].position[0], currTri.v[1].position[1] - currTri.v[0].position[1], currTri.v[1].position[2] - currTri.v[0].position[2]);
          glm::vec3 e2(currTri.v[2].position[0] - currTri.v[0].position[0], currTri.v[2].position[1] - currTri.v[0].position[1], currTri.v[2].position[2] - currTri.v[0].position[2]);
          glm::vec3 e3(currTri.v[2].position[0] - currTri.v[1].position[0], currTri.v[2].position[1] - currTri.v[1].position[1], currTri.v[2].position[2] - currTri.v[1].position[2]);
          glm::vec3 normal1(currTri.v[0].normal[0], currTri.v[0].normal[1], currTri.v[0].normal[2]); // glm::cross(e1, e2);   // from v0
          glm::vec3 normal2(currTri.v[1].normal[0], currTri.v[1].normal[1], currTri.v[1].normal[2]); // glm::cross(e3, -e1);  // from v1
          glm::vec3 normal3(currTri.v[2].normal[0], currTri.v[2].normal[1], currTri.v[2].normal[2]); // glm::cross(-e2, -e3); // from v2
          normal1 = normalize(normal1);
          normal2 = normalize(normal2);
          normal3 = normalize(normal3);

          glm::vec3 intersectfromVertex0(posObjIntersect.x - currTri.v[0].position[0], posObjIntersect.y - currTri.v[0].position[1], posObjIntersect.z - currTri.v[0].position[2]);
          glm::vec3 intersectfromVertex1(posObjIntersect.x - currTri.v[1].position[0], posObjIntersect.y - currTri.v[1].position[1], posObjIntersect.z - currTri.v[1].position[2]);
          float alpha = glm::length(glm::cross(intersectfromVertex1, e3)) / glm::length(glm::cross(e1, e2));
          float beta = glm::length(glm::cross(e2, intersectfromVertex0)) / glm::length(glm::cross(e1, e2));
          float gamma = glm::length(glm::cross(intersectfromVertex0, e1)) / glm::length(glm::cross(e1, e2));
          glm::vec3 normal = normal1 * alpha + beta * normal2 + gamma * normal3;
          normal = normalize(normal);
          // printf("\nnormal: %f, %f, %f", normal.x, normal.y, normal.z);
          glm::vec3 R = (2.0f * glm::dot(objToLight, normal) * normal) - objToLight;
          R = normalize(R);
          //    clamp dot products to 0
          float kdx = float(currTri.v[0].color_diffuse[0] * alpha + beta * currTri.v[1].color_diffuse[0] + gamma * currTri.v[2].color_diffuse[0]);
          float kdy = float(currTri.v[0].color_diffuse[1] * alpha + beta * currTri.v[1].color_diffuse[1] + gamma * currTri.v[2].color_diffuse[1]);
          float kdz = float(currTri.v[0].color_diffuse[2] * alpha + beta * currTri.v[1].color_diffuse[2] + gamma * currTri.v[2].color_diffuse[2]);

          float ksx = float(currTri.v[0].color_specular[0] * alpha + beta * currTri.v[1].color_specular[0] + gamma * currTri.v[2].color_specular[0]);
          float ksy = float(currTri.v[0].color_specular[1] * alpha + beta * currTri.v[1].color_specular[1] + gamma * currTri.v[2].color_specular[1]);
          float ksz = float(currTri.v[0].color_specular[2] * alpha + beta * currTri.v[1].color_specular[2] + gamma * currTri.v[2].color_specular[2]);
          // printf(" alpha: %f, beta: %f, gamma: %f", alpha, beta, gamma);
          float triShine = alpha * currTri.v[0].shininess + beta * currTri.v[1].shininess + gamma * currTri.v[2].shininess;
          float LdotN = glm::dot(objToLight, normal);

          glm::vec3 v = normalize(-curr.direction);
          float RdotV = glm::dot(v, R);

          if (LdotN < 0)
            LdotN = 0.0;
          else if (LdotN > 1.0)
            LdotN = 1.0;
          if (RdotV < 0)
            RdotV = 0.0;
          else if (RdotV > 1.0)
            RdotV = 1.0;

          float RVdotExpo = pow(RdotV, triShine);
          float Ix = float(lights[l].color[0]) * (kdx * LdotN + ksx * RVdotExpo);
          float Iy = float(lights[l].color[1]) * (kdy * LdotN + ksy * RVdotExpo);
          float Iz = float(lights[l].color[2]) * (kdz * LdotN + ksz * RVdotExpo);
          // glm::vec3 c(Ix, Iy, Iz);
          color.x += Ix;
          color.y += Iy;
          color.z += Iz;

          if (isnan(Ix) || isnan(Iy) || isnan(Iz))
          {
            // printf("\nxval: %f and y val: %f and z val: %f", Ix, Iy, Iz);
            // printf("\nkdx: %f and LdotN: %f and ksx: %f and RV: %f", kdx, LdotN, ksx, RVdotExpo);
            printf("\ndot product %f and shine: %f and finalval %f", glm::dot(v, R), triShine, pow(glm::dot(v, R), triShine));
            color.x = 0.5;
            color.y = 0.5;
            color.z = 0.5;
          }
          else if (Ix == 0 || Iy == 0 || Iz == 0)
          {
            // printf("\nxval: %f and y val: %f and z val: %f", Ix, Iy, Iz);
          }
        }
      }
      else
      {
        // printf("\nxval: %f and y val: %f", posObjIntersect.x, posObjIntersect.y);
        // object is in SHADOW!
        color.x = 0;
        color.y = 0;
        color.z = 0;
        // if (intersectObjType == "triangle") // shadow ray of never intersects with triangle
        //  printf("\nshadow ray value: %f, %f, %f", posShadowIntersect.x, posShadowIntersect.y, posShadowIntersect.z);
        // if (intersectObjType == "sphere") // shadow ray of never intersects with triangle
        //   printf("\nshadow ray value: %f, %f, %f", posShadowIntersect.x, posShadowIntersect.y, posShadowIntersect.z);
      }
    }
    // add global ambiant light, clamp to 1.0 for each rbg value
    color.x += float(ambient_light[0]);
    color.y += float(ambient_light[1]);
    color.z += float(ambient_light[2]);
  }
  else
  {
    color.x = 1.0;
    color.y = 1.0;
    color.z = 1.0;
  }
  if (color.x > 1.0)
    color.x = 1.0;
  if (color.y > 1.0)
    color.y = 1.0;
  if (color.z > 1.0)
    color.z = 1.0;
  return color;
  // get normal of intersection point
  // cast shadow rays to each of the light sources
  // throw ambiant light colors at the end
}

float intersect(Ray curr, float lightT)
{
  vector<float> tValuesSphere;
  vector<int> indexesSphere;
  float closeSphereT = MAXFLOAT;
  float closeSphereIndex;

  vector<float> tValuesTri;
  vector<int> indexesTri;
  float closeTriT = MAXFLOAT;
  float closeTriIndex;

  for (int k = 0; k < num_spheres; k++)
  {
    Sphere currSphere = spheres[k];
    float a = dot(curr.direction, curr.direction);
    glm::vec3 originToCenter(curr.origin.x - currSphere.position[0], curr.origin.y - currSphere.position[1], curr.origin.z - currSphere.position[2]);
    float b = 2.0 * dot(originToCenter, curr.direction);
    float c = dot(originToCenter, originToCenter) - (currSphere.radius * currSphere.radius);
    // printf(" direction- a: %f b: %f c: %f \n ", a, b, c);
    // printf(" ray origin - x: %f y: %f z: %f ,  ", originToCenter.x, originToCenter.y, originToCenter.z);
    //   printf("radius %f pos %f %f %f\n", currSphere.radius, currSphere.position[0], currSphere.position[1], currSphere.position[2]);
    float check = b * b - 4 * c * a;

    if ((b * b - 4 * c * a) >= 0) // valid values to get t values
    {
      // printf(" a: %f b: %f c: %f ,  ", a, b, c);
      float discrim1 = (-b + sqrt(b * b - 4 * c * a)) / 2 * a;
      float discrim2 = (-b - sqrt(b * b - 4 * c * a)) / 2 * a;
      // if intersects, do something
      if (discrim1 > 1e-3 || discrim2 > 1e-3)
      {
        float minT = min(discrim1, discrim2);
        tValuesSphere.push_back(minT);
        indexesSphere.push_back(k);
      }
    }
    // printf(" t value: %f  ", tValues.at(0));
  }
  // get closest sphere and its index within Spheres array
  if (tValuesSphere.size() != 0)
  {
    for (int p = 0; p < tValuesSphere.size(); p++)
    {
      if (tValuesSphere.at(p) < closeSphereT)
      {
        closeSphereT = tValuesSphere.at(p);
        closeSphereIndex = indexesSphere.at(p);
      }
    }
  }

  // triangle intersection
  for (int w = 0; w < num_triangles; w++)
  {
    bool checkingShadow = false;
    // if (w == 3)
    //  printf("\nin last triangle");
    Triangle currTri = triangles[w];
    glm::vec3 v0Tov1(currTri.v[1].position[0] - currTri.v[0].position[0], currTri.v[1].position[1] - currTri.v[0].position[1], currTri.v[1].position[2] - currTri.v[0].position[2]);
    glm::vec3 v0Tov2(currTri.v[2].position[0] - currTri.v[0].position[0], currTri.v[2].position[1] - currTri.v[0].position[1], currTri.v[2].position[2] - currTri.v[0].position[2]);
    glm::vec3 v1Tov2(currTri.v[2].position[0] - currTri.v[1].position[0], currTri.v[2].position[1] - currTri.v[1].position[1], currTri.v[2].position[2] - currTri.v[1].position[2]);
    glm::vec3 n = glm::cross(v0Tov1, v0Tov2);
    // v0Tov1 = normalize(v0Tov1);
    // v0Tov2 = normalize(v0Tov2);
    // v1Tov2 = normalize(v1Tov2);
    n = normalize(n);
    glm::vec3 pointOnPlane(currTri.v[0].position[0], currTri.v[0].position[1], currTri.v[0].position[2]);
    // pointOnPlane = normalize(pointOnPlane);
    //  pointOnPlane = normalize(pointOnPlane);

    float dotProductND = glm::dot(n, curr.direction);
    bool notParallel = true;
    if (dotProductND == 0)
    {
      notParallel = false;
    }
    if (notParallel)
    {
      float dCoeff = -glm::dot(n, pointOnPlane);
      // printf("dot : %f , dCoeff: %f", dot(n, curr.origin), dCoeff);
      float t = -(dot(n, curr.origin) + dCoeff) / (dotProductND);
      glm::vec3 intersectionPoint(curr.origin.x + t * curr.direction.x, curr.origin.y + t * curr.direction.y, curr.origin.z + t * curr.direction.z);
      // intersectionPoint = normalize(intersectionPoint);
      glm::vec3 intersectfromVertex0(intersectionPoint.x - currTri.v[0].position[0], intersectionPoint.y - currTri.v[0].position[1], intersectionPoint.z - currTri.v[0].position[2]);
      glm::vec3 intersectfromVertex1(intersectionPoint.x - currTri.v[1].position[0], intersectionPoint.y - currTri.v[1].position[1], intersectionPoint.z - currTri.v[1].position[2]);

      float alpha = glm::length(glm::cross(intersectfromVertex1, v1Tov2)) / glm::length(glm::cross(v0Tov1, v0Tov2));
      float beta = glm::length(glm::cross(v0Tov2, intersectfromVertex0)) / glm::length(glm::cross(v0Tov1, v0Tov2));
      float gamma = glm::length(glm::cross(intersectfromVertex0, v0Tov1)) / glm::length(glm::cross(v0Tov1, v0Tov2));

      bool outside;
      // float smallVal = 1e-5;
      if (alpha >= 0 && alpha < 1.0001f && beta >= 0 && beta < 1.0001f && gamma >= 0 && gamma < 1.0001f && (alpha + beta + gamma) < 1.001f)
      {
        // printf("\nthis is beta: %f and this is alpha %f ", beta, alpha);
        outside = false;
      }
      else
      {
        outside = true;
      }
      if (t > 1e-5 && !outside && t < lightT) // determine if it is intersected !
      {
        // printf("\nthis is the t value: %f ", t);
        tValuesTri.push_back(t);
        indexesTri.push_back(w);
      }
    }
  }
  // get closest t triangle and its index within Triangle array
  if (tValuesTri.size() != 0)
  {
    for (int p = 0; p < tValuesTri.size(); p++)
    {
      if (tValuesTri.at(p) < closeTriT) //
      {
        closeTriT = tValuesTri.at(p);
        closeTriIndex = indexesTri.at(p);
      }
    }
  }

  if (tValuesSphere.size() == 0 && tValuesTri.size() == 0) // return  intersection point if it exists for either spheres or triangles
  {
    return MAXFLOAT; // vector of max float coordinates
  }
  else
  {
    // exclude object where it came from
    // throw away less than small number
    float bestT;
    if (tValuesSphere.size() != 0 && closeTriT > closeSphereT) // sphere is closest
    {
      // rintf(" %f ", closeTriT - closeSphereT);
      bestT = closeSphereT;
      intersectObjIndex = closeSphereIndex;
      intersectObjType = "sphere";
    }
    else
    {
      bestT = closeTriT;
      intersectObjIndex = closeTriIndex;
      intersectObjType = "triangle";
      // t value has to be less than the intersection point and the light
    }
    // glm::vec3 iPoint(curr.origin.x + bestT * curr.direction.x, curr.origin.y + bestT * curr.direction.y, curr.origin.z + bestT * curr.direction.z);
    // glm::vec3 normal(0, 0, 0);
    return bestT; // type of object, index in vertex
  }
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x, y, r, g, b);
  if (mode == MODE_JPEG)
    plot_pixel_jpeg(x, y, r, g, b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if (strcasecmp(expected, found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE *file, const char *check, double p[3])
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check(check, str);
  fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
  printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check("rad:", str);
  fscanf(file, "%lf", r);
  printf("rad: %f\n", *r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file, "%s", s);
  parse_check("shi:", s);
  fscanf(file, "%lf", shi);
  printf("shi: %f\n", *shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv, "r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file, "%i", &number_of_objects);

  printf("number of objects: %i\n", number_of_objects);

  parse_doubles(file, "amb:", ambient_light);

  for (int i = 0; i < number_of_objects; i++)
  {
    fscanf(file, "%s\n", type);
    printf("%s\n", type);
    if (strcasecmp(type, "triangle") == 0)
    {
      printf("found triangle\n");
      for (int j = 0; j < 3; j++)
      {
        parse_doubles(file, "pos:", t.v[j].position);
        parse_doubles(file, "nor:", t.v[j].normal);
        parse_doubles(file, "dif:", t.v[j].color_diffuse);
        parse_doubles(file, "spe:", t.v[j].color_specular);
        parse_shi(file, &t.v[j].shininess);
      }

      if (num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if (strcasecmp(type, "sphere") == 0)
    {
      printf("found sphere\n");

      parse_doubles(file, "pos:", s.position);
      parse_rad(file, &s.radius);
      parse_doubles(file, "dif:", s.color_diffuse);
      parse_doubles(file, "spe:", s.color_specular);
      parse_shi(file, &s.shininess);

      if (num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if (strcasecmp(type, "light") == 0)
    {
      printf("found light\n");
      parse_doubles(file, "pos:", l.position);
      parse_doubles(file, "col:", l.color);

      if (num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n", type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // hack to make it only draw once
  static int once = 0;
  if (!once)
  {
    draw_scene();
    if (mode == MODE_JPEG)
      save_jpg();
  }
  once = 1;
}

int main(int argc, char **argv)
{
  if ((argc < 2) || (argc > 3))
  {
    printf("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if (argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if (argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc, argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(WIDTH, HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
#ifdef __APPLE__
  // This is needed on recent Mac OS X versions to correctly display the window.
  glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
#endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
