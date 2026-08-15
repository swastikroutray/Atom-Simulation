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

float orbitDistance = 15.0f;

struct Engine {
    GLFWwindow* window;
    int WIDTH = 800, HEIGHT = 600;

    Engine() {
        if (!glfwInit()) {
            cerr << "Failed to init GLFW\n";
            exit(EXIT_FAILURE);
        }

        window = glfwCreateWindow(WIDTH, HEIGHT, "2D Atom Sim", nullptr, nullptr);
        if (!window) {
            cerr << "Failed to create window\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);
        
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) {
            cerr << "Note: GLEW failed or not strictly required for fixed pipeline\n";
        }

        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glViewport(0, 0, fbWidth, fbHeight);
    }

    void run() {
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

// Waves
struct WavePoint { 
    vec2 localPos; 
    vec2 dir; 
};

struct Wave {
    float energy;
    float sigma = 40.0f, k = 0.4f, phase = 0.0f, a = 10.0f;
    vector<WavePoint> points;
    vec2 pos, dir;
    vec3 col;

    Wave(float e, vec2 pos, vec2 dir, vec3 col = vec3(0.0f, 1.0f, 1.0f)) 
        : energy(e), pos(pos), dir(dir), col(col) {
        this->dir = normalize(dir);
        // Step increased from 0.1f to 1.0f for smooth 60 FPS on CPU laptops
        for (float x = -sigma; x <= sigma; x += 1.0f) {
            points.push_back({ pos + x * this->dir, this->dir * 200.0f });
        }
    }

    void draw() {
        glColor3f(col.r, col.g, col.b);
        glBegin(GL_LINE_STRIP);
        for (WavePoint& p : points) {
            vec2 perp(-p.dir.y, p.dir.x);
            perp = normalize(perp);
            float y_disp = a * sin(k * length(p.localPos) - phase);
            vec2 drawPos = p.localPos + perp * y_disp;
            glVertex2f(drawPos.x, drawPos.y);
        }
        glEnd();
    }

    bool update(float dt) {
        phase += 30.0f * dt;
        for (WavePoint& p : points) {
            p.localPos += p.dir * dt;
            if (p.localPos.x < -engine.WIDTH / 2.0f || p.localPos.x > engine.WIDTH / 2.0f ||
                p.localPos.y < -engine.HEIGHT / 2.0f || p.localPos.y > engine.HEIGHT / 2.0f) {
                return true; // Out of bounds
            }
        }
        return false;
    }
};

vector<Wave> waves;
vector<Wave> pendingWaves; // Buffer to prevent vector re-allocation during loops

// Particles
struct Particle {
    vec2 pos;
    int charge;
    float angle = 0.0f;
    int n = 1;
    float excitedTimer = 0.0f;

    Particle(vec2 pos, int charge) : pos(pos), charge(charge) {}

    void draw(vec2 centre, int segments = 32) {
        if (charge == -1) {
            glLineWidth(0.4f);
            glBegin(GL_LINE_LOOP);
            glColor3f(0.18f, 0.22f, 0.3f);
            for (int i = 0; i <= segments; i++) {
                float theta = 2.0f * M_PI * i / segments;
                float x = cos(theta) * n * orbitDistance;
                float y = sin(theta) * n * orbitDistance;
                glVertex2f(x + centre.x, y + centre.y);
            }
            glEnd();
        }

        float r = (charge == -1) ? 2.0f : 5.0f;
        if (charge == -1)      glColor3f(0.0f, 0.85f, 1.0f);
        else if (charge == 1)  glColor3f(1.0f, 0.3f, 0.45f);
        else                   glColor3f(0.45f, 0.55f, 0.65f);

        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(pos.x, pos.y);
        for (int i = 0; i <= segments; i++) {
            float theta = 2.0f * M_PI * i / segments;
            float x = cos(theta) * r;
            float y = sin(theta) * r;
            glVertex2f(x + pos.x, y + pos.y);
        }
        glEnd();
    }

    void update(vec2 c) {
        float r = n * orbitDistance;
        angle += 0.05f;
        pos = vec2(cos(angle) * r + c.x, sin(angle) * r + c.y);

        if (excitedTimer <= 0.0f && n > 1) {
            n--;
            excitedTimer = 0.5f; // Reasonable decay interval
            float waveDirX = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
            float waveDirY = ((float)rand() / RAND_MAX) * 2.0f - 1.0f;
            float energyDiff = -13.6f / ((n + 1) * (n + 1)) - (-13.6f / (n * n));
            pendingWaves.emplace_back(energyDiff, pos, normalize(vec2(waveDirX, waveDirY)), vec3(1.0f, 0.85f, 0.2f));
        }
    }
};

struct Atom {
    vec2 pos;
    vec2 v = vec2(0.0f);
    vector<Particle> particles;

    Atom(vec2 p) : pos(p) {
        particles.emplace_back(pos, 1);
        particles.emplace_back(vec2(pos.x - orbitDistance, pos.y), -1);
    }
};

vector<Atom> atoms;

static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button != GLFW_MOUSE_BUTTON_LEFT || action != GLFW_PRESS) return;

    double mx, my;
    glfwGetCursorPos(window, &mx, &my);

    Engine* eng = static_cast<Engine*>(glfwGetWindowUserPointer(window));
    float worldX = (float)mx - eng->WIDTH / 2.0f;
    float worldY = eng->HEIGHT / 2.0f - (float)my;
    vec2 spawnPos(worldX, worldY);

    float energyN1toN2 = -13.6f / 4.0f - (-13.6f);
    for (int i = 0; i < 25; i++) {
        float angle = ((float)rand() / RAND_MAX) * 2.0f * M_PI;
        vec2 dir(cos(angle), sin(angle));
        waves.emplace_back(energyN1toN2, spawnPos, dir);
    }
}

int main() {
    srand(static_cast<unsigned int>(time(nullptr)));

    int num_atoms = 20;
    float radius = 100.0f;
    for (int i = 0; i < num_atoms; i++) {
        float angle = 2.0f * M_PI * i / num_atoms;
        atoms.emplace_back(vec2(cos(angle) * radius, sin(angle) * radius));
    }

    glfwSetWindowUserPointer(engine.window, &engine);
    glfwSetMouseButtonCallback(engine.window, mouseButtonCallback);

    float energyN1toN2 = -13.6f / 4.0f - (-13.6f);
    for (int i = 0; i < 24; i++) {
        waves.emplace_back(energyN1toN2, vec2(200.0f, i * 20 - 200), vec2(-1.0f, 0.0f));
    }

    while (!glfwWindowShouldClose(engine.window)) {
        engine.run();

        // Atom Physics & Drawing
        for (Atom &a : atoms) {
            for (Atom &a2 : atoms) {
                if (&a2 == &a) continue;
                float dist = std::max(length(a.pos - a2.pos), 1.0f);
                vec2 dir = normalize(a.pos - a2.pos);
                a.v += dir / dist * 57.5f;
            }

            // Boundary repulsion
            const float boundary_stiffness = 0.01f;
            const float boundary_threshold = 200.0f;

            float dist_left = a.pos.x + engine.WIDTH / 2.0f;
            if (dist_left < boundary_threshold) a.v.x += (boundary_threshold - dist_left) * boundary_stiffness;

            float dist_right = engine.WIDTH / 2.0f - a.pos.x;
            if (dist_right < boundary_threshold) a.v.x -= (boundary_threshold - dist_right) * boundary_stiffness;

            float dist_top = engine.HEIGHT / 2.0f - a.pos.y;
            if (dist_top < boundary_threshold) a.v.y -= (boundary_threshold - dist_top) * boundary_stiffness;

            float dist_bottom = a.pos.y + engine.HEIGHT / 2.0f;
            if (dist_bottom < boundary_threshold) a.v.y += (boundary_threshold - dist_bottom) * boundary_stiffness;

            a.pos += a.v;
            a.v *= 0.99f;

            for (Particle &p : a.particles) {
                p.draw(a.pos);

                if (p.charge == 1) {
                    p.pos = a.pos;
                } else if (p.charge == -1) {
                    if (p.excitedTimer > 0.0f) {
                        p.excitedTimer -= 0.01f;
                    }
                    p.update(a.pos);

                    // Wave absorption check
                    for (Wave& wave : waves) {
                        if (wave.energy <= 0.0f) continue;

                        float energyforUp = -13.6f / ((p.n + 1) * (p.n + 1)) - (-13.6f / (p.n * p.n));
                        if (fabs(wave.energy - energyforUp) < 0.001f) {
                            for (WavePoint& wp : wave.points) {
                                if (length(p.pos - wp.localPos) < 20.0f) {
                                    wave.energy = 0.0f;
                                    p.n += 1;
                                    p.excitedTimer = 1.0f;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Push newly spawned de-excitation waves safely
        if (!pendingWaves.empty()) {
            waves.insert(waves.end(), pendingWaves.begin(), pendingWaves.end());
            pendingWaves.clear();
        }

        // Draw and Update Waves
        for (auto it = waves.begin(); it != waves.end(); ) {
            if (it->energy <= 0.0f) {
                it = waves.erase(it);
                continue;
            }

            it->draw();
            if (it->update(0.01f)) {
                it = waves.erase(it);
            } else {
                ++it;
            }
        }

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}