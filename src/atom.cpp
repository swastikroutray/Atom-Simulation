#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;

//variables

float orbitDistance = 15.0f;

//engine

struct Wave;
vec2 mouseWorld(0.0f);
struct Engine {

    GLFWwindow* window;
    int WIDTH = 800, HEIGHT = 600;

    Engine () {
        // Init GLFW 
        if (!glfwInit()) {
            cerr << "failed to init glfw, LOL";
            exit(EXIT_FAILURE);
        }

        // Create Window
        window = glfwCreateWindow(WIDTH, HEIGHT, "2D atom sim by swastik", nullptr, nullptr);
        if (!window) {
            cerr << "failed to create window, LOLOLOL";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glfwMakeContextCurrent(window);
        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glViewport(0, 0, fbWidth, fbHeight);
    }
    void run() {
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        // set origin to centre
        double halfWidth = WIDTH / 2.0f, halfHeight = HEIGHT / 2.0f;
        glOrtho(-halfWidth, halfWidth, -halfHeight, halfHeight, -1.0, 1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
};
Engine engine;

//waves

struct WavePoint { vec2 localPos; vec2 dir;  };
struct Wave {
   float energy;
   float sigma = 40.0f, k = 0.4f, phase = 0.0f, a = 10.0f, angleR;
   vector<WavePoint> points;
    vec2 pos, dir;
    vec3 col;
    Wave(float e, vec2 pos, vec2 dir, vec3 col = vec3(0.0f, 1.0f, 1.0f)) : energy(e), pos(pos), dir(dir), col(col) {
       dir = normalize(dir);
       for (float x = -sigma; x <= sigma; x += 0.1f) {
           points.push_back({ pos + x*dir, dir * 200.0f});
       }
       angleR = atan2(dir.y, dir.x);
   }
  
   void draw() {
        glColor3f(col.r, col.g, col.b);
        glBegin(GL_LINE_STRIP);
        for (WavePoint& p : points) {
            // Perpendicular vector for sine displacement
            vec2 perp(-p.dir.y, p.dir.x);
            perp = normalize(perp);


        // Use global phase, not x, for sine
        float y_disp = a * sin(k * length(p.localPos) - phase);


            vec2 drawPos = p.localPos + perp * y_disp;
            glVertex2f(drawPos.x, drawPos.y);
        }
        glEnd();
    }
    bool update(float dt) {
       phase += 30.0f * dt; // continuous phase


       for (WavePoint& p : points) {
           // move along velocity
           p.localPos += p.dir * dt;


            if (p.localPos.x < -engine.WIDTH/2.0f || p.localPos.x > engine.WIDTH/2.0f || p.localPos.y < -engine.HEIGHT/2.0f || p.localPos.y > engine.HEIGHT/2.0f) {
                return true;
            }
       }
       return false;
   }
};
vector<Wave> waves { };

//particles
