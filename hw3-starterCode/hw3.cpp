/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Rachel Bakke
 * *************************
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

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);
glm::vec3 colorRay(Ray curr);
float intersect(Ray curr, float lightT);

glm::vec3 maxVal(MAXFLOAT, MAXFLOAT, MAXFLOAT);

void draw_scene()
{
  int myHeight = 480;
  int myWidth = 640;

  float aVal = (myWidth * 1.0f) / (myHeight * 1.0f);
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
  float wStepSize = (width3D * 1.0f) / (1.0f * myWidth);
  float hStepSize = (height3D * 1.0f) / (myHeight * 1.0f);

  // PLOTTING PIXELS - get color of the pixel from the ray
  for (unsigned int x = 0; x < myWidth; x++)
  {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for (unsigned int y = 0; y < myHeight; y++)
    {
      glm::vec3 finalColor(0, 0, 0);
      for (int j = 0; j < 2; j++)
      {
        for (int i = 0; i < 2; i++)
        {
          float xVal = (x * wStepSize) + botLeftX;
          float yVal = (y * hStepSize) + botLeftY;
          glm::vec3 direction(xVal + i * 0.001, yVal + j * 0.001, -1.0);
          direction = normalize(direction);
          glm::vec3 origin1(0, 0, 0);
          Ray currRay = Ray(origin1, direction);

          glm::vec3 subpixelColor = colorRay(currRay);
          finalColor.x += subpixelColor.x;
          finalColor.y += subpixelColor.y;
          finalColor.z += subpixelColor.z;
        }
      }
      plot_pixel(x, y, finalColor.x / 4.0f * 255 * 0.7, finalColor.y / 4.0f * 255 * 0.7, finalColor.z / 4.0f * 255 * 0.7);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n");
  fflush(stdout);
}

glm::vec3 colorRay(Ray curr)
{
  glm::vec3 color(0, 0, 0); // set as black
  float currT = intersect(curr, MAXFLOAT);
  glm::vec3 posObjIntersect(curr.origin.x + currT * curr.direction.x, curr.origin.y + currT * curr.direction.y, curr.origin.z + currT * curr.direction.z);
  bool checkingTri = false;
  if (currT != MAXFLOAT) // object intersects
  {
    // check shadow rays for each light source and see if they are hit
    for (int l = 0; l < num_lights; l++) //
    {
      // create more lights, loop through those lights,
      Ray currShadow = Ray();
      currShadow.setOrigin(posObjIntersect);
      glm::vec3 objToLight(lights[l].position[0] - posObjIntersect.x, lights[l].position[1] - posObjIntersect.y, lights[l].position[2] - posObjIntersect.z);
      float distanceToLight = glm::length(objToLight);
      objToLight = normalize(objToLight);

      glm::vec3 newShadowOrigin(posObjIntersect.x + 0.00001f * objToLight.x, posObjIntersect.y + 0.00001f * objToLight.y, posObjIntersect.z + 0.00001f * objToLight.z);
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

      if (shadowT == MAXFLOAT && shadowT > 1e-3) // NOT IN SHADOW
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

          float kdx = float(currSphere.color_diffuse[0]);
          float kdy = float(currSphere.color_diffuse[1]);
          float kdz = float(currSphere.color_diffuse[2]);

          float ksx = float(currSphere.color_specular[0]);
          float ksy = float(currSphere.color_specular[1]);
          float ksz = float(currSphere.color_specular[2]);

          glm::vec3 v = normalize(-curr.direction);
          //     clamp dot products to 0
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
          glm::vec3 R = (2.0f * glm::dot(objToLight, normal) * normal) - objToLight;
          R = normalize(R);

          float kdx = float(currTri.v[0].color_diffuse[0] * alpha + beta * currTri.v[1].color_diffuse[0] + gamma * currTri.v[2].color_diffuse[0]);
          float kdy = float(currTri.v[0].color_diffuse[1] * alpha + beta * currTri.v[1].color_diffuse[1] + gamma * currTri.v[2].color_diffuse[1]);
          float kdz = float(currTri.v[0].color_diffuse[2] * alpha + beta * currTri.v[1].color_diffuse[2] + gamma * currTri.v[2].color_diffuse[2]);

          float ksx = float(currTri.v[0].color_specular[0] * alpha + beta * currTri.v[1].color_specular[0] + gamma * currTri.v[2].color_specular[0]);
          float ksy = float(currTri.v[0].color_specular[1] * alpha + beta * currTri.v[1].color_specular[1] + gamma * currTri.v[2].color_specular[1]);
          float ksz = float(currTri.v[0].color_specular[2] * alpha + beta * currTri.v[1].color_specular[2] + gamma * currTri.v[2].color_specular[2]);
          float triShine = alpha * currTri.v[0].shininess + beta * currTri.v[1].shininess + gamma * currTri.v[2].shininess;
          float LdotN = glm::dot(objToLight, normal);

          glm::vec3 v = normalize(-curr.direction);
          float RdotV = glm::dot(v, R);
          //    clamp dot products to 0
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
          // printf("\nlight num: %d colorx: %f, colory: %f, colorz: %f", l, color.x * 255, color.y * 255, color.z);
        }
      }
      else
      {
        // printf("\nshadow + light num: %d colorx: %f, colory: %f, colorz: %f", l, color.x * 255, color.y * 255, color.z);
        // object is in SHADOW!
        // color.x = 0; // color.x * 0.1;
        // color.y = 0; // color.y * 0.1;
        // color.z = 0; // color.z * 0.1;
      }
    } // end of singular light routine
    // add global ambiant light, clamp to 1.0 for each rbg value
  } // end of intersection point
  else
  {
    // blue
    /** color.x = 0.537254902;
     color.y = 0.8117647059;
     color.z = 0.9411764706;*/
    color.x = 1.0;
    color.y = 1.0;
    color.z = 1.0;
  }
  // add ambient light

  color.x += float(ambient_light[0]);
  color.y += float(ambient_light[1]);
  color.z += float(ambient_light[2]);
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
    float check = b * b - 4 * c * a;

    if ((b * b - 4 * c * a) >= 0) // valid values to get t values
    {
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
    Triangle currTri = triangles[w];
    glm::vec3 v0Tov1(currTri.v[1].position[0] - currTri.v[0].position[0], currTri.v[1].position[1] - currTri.v[0].position[1], currTri.v[1].position[2] - currTri.v[0].position[2]);
    glm::vec3 v0Tov2(currTri.v[2].position[0] - currTri.v[0].position[0], currTri.v[2].position[1] - currTri.v[0].position[1], currTri.v[2].position[2] - currTri.v[0].position[2]);
    glm::vec3 v1Tov2(currTri.v[2].position[0] - currTri.v[1].position[0], currTri.v[2].position[1] - currTri.v[1].position[1], currTri.v[2].position[2] - currTri.v[1].position[2]);
    glm::vec3 n = glm::cross(v0Tov1, v0Tov2);
    n = normalize(n);
    glm::vec3 pointOnPlane(currTri.v[0].position[0], currTri.v[0].position[1], currTri.v[0].position[2]);

    float dotProductND = glm::dot(n, curr.direction);
    bool notParallel = true;
    if (dotProductND >= .00001f && dotProductND <= -.00001f)
    {
      notParallel = false;
    }
    if (notParallel)
    {
      float dCoeff = -glm::dot(n, pointOnPlane);
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
      if ((alpha + beta + gamma) < 1.00001f)
      {
        outside = false;
      }
      else
      {
        outside = true;
      }
      if (t > 1e-3 && !outside && t < lightT) // determine if it is intersected !
      {
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
    float bestT;
    if (tValuesSphere.size() != 0 && closeTriT > closeSphereT) // sphere is closest
    {
      bestT = closeSphereT;
      intersectObjIndex = closeSphereIndex;
      intersectObjType = "sphere";
    }
    else
    {
      bestT = closeTriT;
      intersectObjIndex = closeTriIndex;
      intersectObjType = "triangle";
    }

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
  /** int count = 0;
   while (count < 5)
   {

     draw_scene();
     for (int i = 0; i < num_spheres; i++)
     {
       spheres[i].position[1] + 5;
     }
     if (mode == MODE_JPEG)
       save_jpg();
     count++;
     for (int i = 0; i < num_spheres; i++)
  {
    if (i % 2 == 0)
      spheres[i].position[1] -= 7;
    else
    {
      spheres[i].position[1] -= 9;
    }
  }
  for (int i = 0; i < num_spheres; i++)
  {
    if (i % 2 == 0)
      spheres[i].position[1] -= 5;
    else
    {
      spheres[i].position[1] -= 8;
    }
  }
   for (int i = 0; i < num_spheres; i++)
  {
    if (i % 2 == 0)
    {
      spheres[i].position[1] += 2;
      spheres[i].position[1] -= 3;
    }
    else
    {
      spheres[i].position[1] += 1;
      spheres[i].position[1] -= 1;
    }
  }
    for (int i = 0; i < num_spheres; i++)
  {
    if (i % 2 == 0)
    {
      spheres[i].position[1] += 15;
      spheres[i].position[1] -= 2;
    }
    else
    {
      spheres[i].position[1] += 11;
      spheres[i].position[1] += 3;
    }
  }
   } */

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
