/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Rachel Bakke
 * *************************
 * ./hw3 ./spheres.scene output.jpg
 *  ./hw3 ./test1.scene output.jpg
 * ./hw3 ./test2.scene output.jpg
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
glm::vec3 intersectObjPos(0, 0, 0);

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);
glm::vec3 colorRay(Ray curr);
glm::vec3 intersect(Ray curr);

glm::vec3 maxVal(MAXFLOAT, MAXFLOAT, MAXFLOAT);

// MODIFY THIS FUNCTION
void draw_scene()
{

  // get screen setup
  Ray rays[WIDTH][HEIGHT];
  printf("START");
  for (int i = 0; i < WIDTH; i++)
  {
    for (int j = 0; j < HEIGHT; j++)
    {
      // STEP 1! generate rays
      float aVal = (WIDTH * 1.0f) / (HEIGHT * 1.0f);
      float xVal = (2.0f * (i + 0.5f) / (WIDTH * 1.0f) - 1.0f) * aVal * tan(glm::radians(fov) / 2.0f);
      float yVal = (1.0f - 2.0f * (j + 0.5f) / (HEIGHT * 1.0f) * tan(glm::radians(fov) / 2.0f));
      glm::vec3 direction(xVal, -yVal, -1); // todo
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
      // glm::vec3 currColor = colorRay(rays[x][y]);
      // plot_pixel(x, y, currColor.x, currColor.y, currColor.z); //FINAL COLOR
      glm::vec3 intersectPoint = intersect(rays[x][y]);
      // printf(" %f, %f, %f", intersectPoint.x, intersectPoint.y, intersectPoint.z);
      if (intersectPoint != maxVal) // yes intersects
      {
        plot_pixel(x, y, 255, 192, 203); // COLOR INTERSECT
      }
      else
      {
        // printf("NOINTERSECTION");
        plot_pixel(x, y, 0.0f, 0.0f, 0.0f);
      }
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
  glm::vec3 posIntersect = intersect(curr);
  // printf(" TEST %f , %f , %f", posIntersect.x, posIntersect.y, posIntersect.z);
  if (posIntersect != maxVal) // object intersects
  {
    // check shadow rays for each light source and see if they are hit
    for (int l = 0; l < num_lights; l++)
    {
      Ray currShadow = Ray();
      currShadow.setOrigin(posIntersect);
      glm::vec3 directionLight(lights[l].position[0] - posIntersect.x, lights[l].position[1] - posIntersect.y, lights[l].position[2] - posIntersect.z);
      directionLight = normalize(directionLight);
      currShadow.setDirection(directionLight);
      glm::vec3 posShadowIntersect = intersect(currShadow);
      if (posShadowIntersect != maxVal) // make sure there was an object intersected
      {
        // phong model solve with normals
        // get normal of the intersection and values of shininess and everything
        if (intersectObjType == "sphere")
        {
          Sphere currSphere = spheres[intersectObjIndex];
          printf("sphere check %d", intersectObjIndex);
          // need to know shape type and index to access it here and get its normal
          float Ix = lights[l].color[0] * directionLight.x;
          float Iy = lights[l].color[1];
          float Iz = lights[l].color[2];
        }
        else
        {
          Triangle currTriangle = triangles[intersectObjIndex];
        }
      }
    }
    // add global ambiant light, clamp to 1.0 for each rbg value
  }
  return color;
  // get normal of intersection point
  // cast shadow rays to each of the light sources
  // throw ambiant light colors at the end
}

glm::vec3 intersect(Ray curr)
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
    // printf(" direction- x: %f y: %f z: %f \n ", curr.direction.x, curr.direction.y, curr.direction.z);
    //   printf(" ray origin - x: %f y: %f z: %f ,  ", currSphere.position[0], currSphere.position[1], currSphere.position[2]);
    //  printf("radius %f pos %f %f %f\n", currSphere.radius, currSphere.position[0], currSphere.position[1], currSphere.position[2]);
    float check = b * b - 4 * c * a;

    if ((b * b - 4 * c * a) >= 0) // valid values to get t values
    {
      // printf(" a: %f b: %f c: %f ,  ", a, b, c);
      float discrim1 = (-b + sqrt(b * b - 4 * c * a)) / 2 * a;
      float discrim2 = (-b - sqrt(b * b - 4 * c * a)) / 2 * a;
      // if intersects, do something
      if (discrim1 >= 0.0001 || discrim2 >= 0.0001)
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
    Triangle currTri = triangles[w];
    currTri.v[0];
    currTri.v[1];
    currTri.v[2];
    glm::vec3 e1(currTri.v[1].position[0] - currTri.v[0].position[0], currTri.v[1].position[1] - currTri.v[0].position[1], currTri.v[1].position[2] - currTri.v[0].position[2]);
    glm::vec3 e2(currTri.v[2].position[0] - currTri.v[0].position[0], currTri.v[2].position[1] - currTri.v[0].position[1], currTri.v[2].position[2] - currTri.v[0].position[2]);
    glm::vec3 n = glm::cross(e1, e2);
    glm::vec3 q = curr.direction;
    glm::vec3 p(curr.origin.x - currTri.v[0].position[0], curr.origin.y - currTri.v[0].position[1], curr.origin.z - currTri.v[0].position[2]);

    float dot = glm::dot(n, q);
    bool notParallel = true;
    if (dot == 0)
    {
      printf("parallel");
      notParallel = true;
    }
    if (notParallel)
    {

      float d = -glm::dot(n, p);
      float t = d / dot;
      glm::vec3 intersectionPoint = curr.origin + t * q;
      glm::vec3 intersectfromPoint1(intersectionPoint.x - currTri.v[0].position[0], intersectionPoint.y - currTri.v[0].position[1], intersectionPoint.z - currTri.v[0].position[2]);

      float u = glm::dot(glm::cross(e2, intersectfromPoint1), n) / glm::dot(n, n);
      float v = glm::dot(glm::cross(intersectfromPoint1, e1), n) / glm::dot(n, n);
      bool notOutside = true;
      if (u < 0 || v < 0 || u + v > 1)
      {
        notOutside = false;
      }

      if (t >= 0.001 && notOutside) // determine if it is intersected
      {
        // printf("tvalue: %f ", t);
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
      if (tValuesTri.at(p) < closeTriT)
      {
        closeTriT = tValuesTri.at(p);
        closeTriIndex = indexesTri.at(p);
      }
    }
  }

  if ((tValuesSphere.size() != 0 || tValuesTri.size() != 0)) // return  intersection point if it exists for either spheres or triangles
  {
    float bestT;
    // printf("TRIANGLE: %f ", closeTriT);
    // if (closeSphereT < closeTriT)
    //  printf("SPHERE: %f ", closeSphereT);
    if (closeSphereT < closeTriT || closeTriT > 2.0f) // sphere is closest
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
    glm::vec3 iPoint(curr.origin.x + bestT * curr.direction.x, curr.origin.y + bestT * curr.direction.y, curr.origin.z + bestT * curr.direction.z);
    // glm::vec3 normal(0, 0, 0);
    return iPoint; // type of object, index in vertex
  }
  else
  {
    return maxVal; // vector of max float coordinates
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
