#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace glm;
using namespace std;

// Base Bohr radius scale in pixels
const float a0 = 25.0f;

struct Engine {
    GLFWwindow* window;
    int WIDTH = 800, HEIGHT = 600;

    Engine() {
        if (!glfwInit()) {
            cerr << "Failed to initialize GLFW\n";
            exit(EXIT_FAILURE);
        }

        window = glfwCreateWindow(WIDTH, HEIGHT, "2D Bohr - de Broglie Wave Sim", nullptr, nullptr);
        if (!window) {
            cerr << "Failed to create GLFW window\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);
        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glViewport(0, 0, fbWidth, fbHeight);

        // Smooth antialiased lines
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_LINE_SMOOTH);
    }

    void run() {
        glClearColor(0.05f, 0.06f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double halfWidth = WIDTH / 2.0;
        double halfHeight = HEIGHT / 2.0;
        glOrtho(-halfWidth, halfWidth, -halfHeight, halfHeight, -1.0, 1.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
};

Engine engine;

//  Particle 
struct Particle {
    vec2 pos;
    int charge;
    float angle = 0.0f;
    float energy = -13.6f; // Ground state n = 1

    Particle(vec2 pos, int charge) : pos(pos), charge(charge) {}

    // Calculate effective principal quantum number n from energy
    float get_n() const {
        if (energy >= -0.01f) return 5.0f; // Ionized / unbound limit
        return sqrt(-13.6f / energy);
    }

    void draw(vec2 centre) {
        // Draw standing de Broglie wave for bound electron
        if (charge == -1) {
            float n = get_n();
            float baseRadius = a0 * (n * n);
            float amplitude = 6.0f * n;
            int segments = 600;

            // Orbit guide ring
            glLineWidth(0.8f);
            glColor4f(0.25f, 0.30f, 0.38f, 0.4f);
            glBegin(GL_LINE_LOOP);
            for (int i = 0; i <= segments; i++) {
                float theta = 2.0f * M_PI * i / segments;
                glVertex2f(cos(theta) * baseRadius + centre.x, sin(theta) * baseRadius + centre.y);
            }
            glEnd();

            // De Broglie standing matter wave
            glLineWidth(1.8f);
            if (fabs(n - round(n)) < 0.08f) {
                glColor3f(0.1f, 0.95f, 0.7f); // Resonant standing wave (integer n)
            } else {
                glColor3f(0.9f, 0.3f, 0.3f); // Non-resonant / destructive interference
            }

            glBegin(GL_LINE_LOOP);
            for (int i = 0; i <= segments; i++) {
                float theta = 2.0f * M_PI * i / segments;
                float r = baseRadius + amplitude * sin(n * theta);
                float x = cos(theta) * r;
                float y = sin(theta) * r;
                glVertex2f(x + centre.x, y + centre.y);
            }
            glEnd();
        }

        // Draw particle sphere
        float r = (charge == -1) ? 6.0f : 12.0f;
        int circle_segments = 32;

        if (charge == -1)      glColor3f(0.0f, 0.85f, 1.0f); // Electron
        else if (charge == 1)  glColor3f(1.0f, 0.25f, 0.35f); // Proton
        else                   glColor3f(0.6f, 0.6f, 0.6f);

        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(pos.x, pos.y);
        for (int i = 0; i <= circle_segments; i++) {
            float theta = 2.0f * M_PI * i / circle_segments;
            glVertex2f(pos.x + cos(theta) * r, pos.y + sin(theta) * r);
        }
        glEnd();
    }

    void update(vec2 c) {
        float n = get_n();
        float baseRadius = a0 * (n * n);
        float amplitude = 6.0f * n;
        float r = baseRadius + amplitude * sin(n * angle);

        // Angular velocity decreases for higher orbits (Keplerian / Bohr speed)
        float angular_speed = 0.04f / (n * 0.75f);
        angle += angular_speed;

        pos = vec2(cos(angle) * r + c.x, sin(angle) * r + c.y);
    }
};

//  Atom 
struct Atom {
    vec2 pos;
    vector<Particle> particles;

    Atom(vec2 p) : pos(p) {
        particles.emplace_back(pos, 1);  // Nucleus
        particles.emplace_back(vec2(pos.x - a0, pos.y), -1); // Electron
    }
};

vector<Atom> atoms = { Atom(vec2(0.0f, 0.0f)) };

// Adjust electron energy
void adjustEnergy(float delta) {
    for (Atom &a : atoms) {
        for (Particle &p : a.particles) {
            if (p.charge == -1) {
                p.energy = std::clamp(p.energy + delta, -13.6f, -0.2f);
                cout << "Energy: " << p.energy << " eV | n = " << p.get_n() << endl;
            }
        }
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action == GLFW_PRESS || action == GLFW_REPEAT) {
        if (key == GLFW_KEY_W) adjustEnergy(0.02f);
        if (key == GLFW_KEY_S) adjustEnergy(-0.02f);
        if (key == GLFW_KEY_E) adjustEnergy(0.2f);
        if (key == GLFW_KEY_D) adjustEnergy(-0.2f);
        if (key == GLFW_KEY_R) adjustEnergy(1.0f);
        if (key == GLFW_KEY_F) adjustEnergy(-1.0f);

        // Quick presets for pure quantum states (n = 1, 2, 3, 4)
        if (key == GLFW_KEY_1) adjustEnergy(-13.60f - atoms[0].particles[1].energy);
        if (key == GLFW_KEY_2) adjustEnergy((-13.6f / 4.0f) - atoms[0].particles[1].energy);
        if (key == GLFW_KEY_3) adjustEnergy((-13.6f / 9.0f) - atoms[0].particles[1].energy);
        if (key == GLFW_KEY_4) adjustEnergy((-13.6f / 16.0f) - atoms[0].particles[1].energy);
    }
}

int main() {
    glfwSetKeyCallback(engine.window, key_callback);

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        for (Atom &a : atoms) {
            for (Particle &p : a.particles) {
                if (p.charge == 1) {
                    p.pos = a.pos;
                } else if (p.charge == -1) {
                    p.update(a.pos);
                }
                p.draw(a.pos);
            }
        }

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}