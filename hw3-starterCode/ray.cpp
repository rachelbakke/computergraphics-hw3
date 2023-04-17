#include <glm/glm.hpp>
using namespace std;

class Ray
{                        // The class
public:                  // Access specifier
    glm::vec3 direction; // Attribute (int variable)
    glm::vec3 origin;    // Attribute (string variable)
    Ray()
    {
        glm::vec3 origin2((0, 0, 0));
        glm::vec3 origin3((0, 0, 0));
        setOrigin(origin2);
        setDirection(origin3);
    }
    Ray(glm::vec3 o, glm::vec3 d)
    {
        setOrigin(o);
        setDirection(d);
    }
    void setOrigin(glm::vec3 input)
    {
        origin = input;
    }
    void setDirection(glm::vec3 input)
    {
        direction = input;
    }
};