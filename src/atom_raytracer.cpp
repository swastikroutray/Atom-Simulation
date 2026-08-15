#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace glm;
using namespace std;

// Physics Constants
const float a0 = 1.0f;
const double hbar = 1.0;
const double m_e = 1.0;

// Simulation State
int N = 50000; // Balanced for CPU/iGPU real-time simulation
int n = 3;
int l = 1;
int m = 1;

struct Particle {
    vec3 pos;
    vec3 vel = vec3(0.0f);
    vec4 color;
    Particle(vec3 p, vec4 c) : pos(p), color(c) {}
};

vector<Particle> particles;
random_device rd;
mt19937 gen(rd());

// Quantum Mechanical Sampling Functions (Dynamic CDFs)
double sampleR(int n, int l, mt19937& rng) {
    const int SAMPLES = 2048;
    const double rMax = 12.0 * n * n * a0;
    vector<double> cdf(SAMPLES);
    double dr = rMax / (SAMPLES - 1);
    double sum = 0.0;

    for (int i = 0; i < SAMPLES; ++i) {
        double r = i * dr;
        double rho = 2.0 * r / (n * a0);

        int k = n - l - 1;
        int alpha = 2 * l + 1;

        double L = 1.0, Lm1 = 1.0 + alpha - rho;
        if (k == 1) L = Lm1;
        else if (k > 1) {
            double Lm2 = 1.0;
            for (int j = 2; j <= k; ++j) {
                L = ((2 * j - 1 + alpha - rho) * Lm1 - (j - 1 + alpha) * Lm2) / j;
                Lm2 = Lm1;
                Lm1 = L;
            }
        }

        double norm = pow(2.0 / (n * a0), 3) * tgamma(n - l) / (2.0 * n * tgamma(n + l + 1));
        double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l) * L;
        double pdf = r * r * R * R;
        sum += pdf;
        cdf[i] = sum;
    }

    for (double& v : cdf) v /= (sum > 0.0 ? sum : 1.0);

    uniform_real_distribution<double> dis(0.0, 1.0);
    double u = dis(rng);
    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * dr;
}

double sampleTheta(int l, int m_val, mt19937& rng) {
    const int SAMPLES = 1024;
    int abs_m = std::abs(m_val);
    vector<double> cdf(SAMPLES);
    double dtheta = M_PI / (SAMPLES - 1);
    double sum = 0.0;

    for (int i = 0; i < SAMPLES; ++i) {
        double theta = i * dtheta;
        double x = cos(theta);

        double Pmm = 1.0;
        if (abs_m > 0) {
            double somx2 = sqrt(std::max(0.0, (1.0 - x) * (1.0 + x)));
            double fact = 1.0;
            for (int j = 1; j <= abs_m; ++j) {
                Pmm *= -fact * somx2;
                fact += 2.0;
            }
        }

        double Plm = Pmm;
        if (l == abs_m) {
            Plm = Pmm;
        } else {
            double Pm1m = x * (2 * abs_m + 1) * Pmm;
            if (l == abs_m + 1) {
                Plm = Pm1m;
            } else {
                for (int ll = abs_m + 2; ll <= l; ++ll) {
                    double Pll = ((2 * ll - 1) * x * Pm1m - (ll + abs_m - 1) * Pmm) / (ll - abs_m);
                    Pmm = Pm1m;
                    Pm1m = Pll;
                }
                Plm = Pm1m;
            }
        }
        double pdf = sin(theta) * Plm * Plm;
        sum += pdf;
        cdf[i] = sum;
    }

    for (double& v : cdf) v /= (sum > 0.0 ? sum : 1.0);

    uniform_real_distribution<double> dis(0.0, 1.0);
    double u = dis(rng);
    int idx = lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    return idx * dtheta;
}

vec4 heatmap_fire(float value) {
    value = std::clamp(value, 0.0f, 1.0f);
    const int num_stops = 6;
    vec4 colors[num_stops] = {
        {0.05f, 0.0f, 0.1f, 0.3f},  // Deep Background Glow
        {0.35f, 0.0f, 0.65f, 0.6f}, // Violet
        {0.85f, 0.1f, 0.0f, 0.8f},  // Deep Red
        {1.0f, 0.55f, 0.0f, 0.9f},  // Orange
        {1.0f, 0.95f, 0.2f, 1.0f},  // Bright Yellow
        {1.0f, 1.0f, 1.0f, 1.0f}    // Core White
    };
    float scaled_v = value * (num_stops - 1);
    int i = static_cast<int>(scaled_v);
    int next_i = std::min(i + 1, num_stops - 1);
    float t = scaled_v - i;
    return mix(colors[i], colors[next_i], t);
}

vec4 calculateOrbitalColor(double r, double theta, int n_val, int l_val, int m_val) {
    double rho = 2.0 * r / (n_val * a0);
    int k = n_val - l_val - 1;
    int alpha = 2 * l_val + 1;

    double L = 1.0, Lm1 = 1.0 + alpha - rho;
    if (k == 1) L = Lm1;
    else if (k > 1) {
        double Lm2 = 1.0;
        for (int j = 2; j <= k; ++j) {
            L = ((2 * j - 1 + alpha - rho) * Lm1 - (j - 1 + alpha) * Lm2) / j;
            Lm2 = Lm1;
            Lm1 = L;
        }
    }
    double norm = pow(2.0 / (n_val * a0), 3) * tgamma(n_val - l_val) / (2.0 * n_val * tgamma(n_val + l_val + 1));
    double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, l_val) * L;
    double radial = R * R;

    double x = cos(theta);
    int abs_m = std::abs(m_val);
    double Pmm = 1.0;
    if (abs_m > 0) {
        double somx2 = sqrt(std::max(0.0, (1.0 - x) * (1.0 + x)));
        double fact = 1.0;
        for (int j = 1; j <= abs_m; ++j) {
            Pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    double Plm = Pmm;
    if (l_val > abs_m) {
        double Pm1m = x * (2 * abs_m + 1) * Pmm;
        if (l_val == abs_m + 1) Plm = Pm1m;
        else {
            for (int ll = abs_m + 2; ll <= l_val; ++ll) {
                double Pll = ((2 * ll - 1) * x * Pm1m - (ll + abs_m - 1) * Pmm) / (ll - abs_m);
                Pmm = Pm1m;
                Pm1m = Pll;
            }
            Plm = Pm1m;
        }
    }
    double intensity = radial * (Plm * Plm);
    float t = static_cast<float>((log10(intensity + 1e-12) + 8.0) / 8.0);
    return heatmap_fire(t);
}

vec3 sphericalToCartesian(float r, float theta, float phi) {
    return vec3(r * sin(theta) * cos(phi), r * cos(theta), r * sin(theta) * sin(phi));
}

vec3 calculateProbabilityFlow(const Particle& p, int m_val) {
    double r = length(p.pos);
    if (r < 1e-5) return vec3(0.0f);
    double theta = acos(std::clamp(p.pos.y / (float)r, -1.0f, 1.0f));
    double phi = atan2(p.pos.z, p.pos.x);

    double sinTheta = sin(theta);
    if (abs(sinTheta) < 1e-3) sinTheta = (sinTheta >= 0 ? 1e-3 : -1e-3);
    double v_mag = (hbar * m_val) / (m_e * r * sinTheta);

    return vec3(-v_mag * sin(phi), 0.0, v_mag * cos(phi));
}

void generateParticles(int count) {
    particles.clear();
    particles.reserve(count);
    uniform_real_distribution<float> disPhi(0.0f, 2.0f * (float)M_PI);

    for (int i = 0; i < count; ++i) {
        double r = sampleR(n, l, gen);
        double theta = sampleTheta(l, m, gen);
        double phi = disPhi(gen);

        vec3 pos = sphericalToCartesian((float)r, (float)theta, (float)phi);
        vec4 col = calculateOrbitalColor(r, theta, n, l, m);
        particles.emplace_back(pos, col);
    }
}
//Interactive Camera
struct Camera {
    vec3 target = vec3(0.0f);
    float radius = 40.0f;
    float azimuth = 0.8f;
    float elevation = 1.2f;
    float orbitSpeed = 0.006f;
    bool dragging = false;
    double lastX = 0.0, lastY = 0.0;

    vec3 position() const {
        float el = std::clamp(elevation, 0.05f, (float)M_PI - 0.05f);
        return vec3(
            radius * sin(el) * cos(azimuth),
            radius * cos(el),
            radius * sin(el) * sin(azimuth)
        );
    }

    void processMouseMove(double x, double y) {
        if (dragging) {
            azimuth += (float)(x - lastX) * orbitSpeed;
            elevation -= (float)(y - lastY) * orbitSpeed;
            elevation = std::clamp(elevation, 0.05f, (float)M_PI - 0.05f);
        }
        lastX = x;
        lastY = y;
    }

    void processMouseButton(int button, int action, GLFWwindow* win) {
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            if (action == GLFW_PRESS) {
                dragging = true;
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if (action == GLFW_RELEASE) {
                dragging = false;
            }
        }
    }

    void processScroll(double yoffset) {
        radius -= (float)yoffset * 2.5f;
        if (radius < 2.0f) radius = 2.0f;
    }
};
Camera camera;

//OpenGl engine and shaders
struct Engine {
    GLFWwindow* window;
    int WIDTH = 1024, HEIGHT = 768;
    GLuint shaderProgram;
    GLuint VAO, VBO;

    Engine() {
        if (!glfwInit()) exit(EXIT_FAILURE);

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

        window = glfwCreateWindow(WIDTH, HEIGHT, "Hydrogen Orbital Wavefunction Visualizer", nullptr, nullptr);
        if (!window) { glfwTerminate(); exit(EXIT_FAILURE); }

        glfwMakeContextCurrent(window);
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK) { glfwTerminate(); exit(EXIT_FAILURE); }

        glEnable(GL_PROGRAM_POINT_SIZE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE); // Additive blending for volumetric glow
        glDisable(GL_DEPTH_TEST);

        initShaders();
        initBuffers();
    }

    void initBuffers() {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, N * sizeof(Particle), nullptr, GL_DYNAMIC_DRAW);

        // Position attribute
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)offsetof(Particle, pos));

        // Color attribute
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)offsetof(Particle, color));

        glBindVertexArray(0);
    }

    void initShaders() {
        const char* vShader = R"(
            #version 330 core
            layout (location = 0) in vec3 aPos;
            layout (location = 1) in vec4 aColor;

            uniform mat4 uMVP;
            uniform float uPointSize;
            out vec4 vColor;

            void main() {
                gl_Position = uMVP * vec4(aPos, 1.0);
                // Perspective point size attenuation
                gl_PointSize = clamp(uPointSize * (60.0 / gl_Position.w), 1.5, 30.0);
                vColor = aColor;
            }
        )";

        const char* fShader = R"(
            #version 330 core
            in vec4 vColor;
            out vec4 FragColor;

            void main() {
                // Soft circular Gaussian point sprite
                vec2 coord = gl_PointCoord - vec2(0.5);
                float distSq = dot(coord, coord);
                if (distSq > 0.25) discard;

                float alpha = exp(-distSq * 10.0);
                FragColor = vec4(vColor.rgb, vColor.a * alpha * 0.35);
            }
        )";

        GLuint vs = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vs, 1, &vShader, nullptr);
        glCompileShader(vs);

        GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fs, 1, &fShader, nullptr);
        glCompileShader(fs);

        shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vs);
        glAttachShader(shaderProgram, fs);
        glLinkProgram(shaderProgram);

        glDeleteShader(vs);
        glDeleteShader(fs);
    }

    void setupCallbacks() {
        glfwSetWindowUserPointer(window, &camera);
        glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int btn, int act, int) {
            ((Camera*)glfwGetWindowUserPointer(win))->processMouseButton(btn, act, win);
        });
        glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
            ((Camera*)glfwGetWindowUserPointer(win))->processMouseMove(x, y);
        });
        glfwSetScrollCallback(window, [](GLFWwindow* win, double, double yoff) {
            ((Camera*)glfwGetWindowUserPointer(win))->processScroll(yoff);
        });
        glfwSetKeyCallback(window, [](GLFWwindow*, int key, int, int act, int) {
            if (act != GLFW_PRESS && act != GLFW_REPEAT) return;
            bool changed = false;

            if (key == GLFW_KEY_W) { n++; changed = true; }
            else if (key == GLFW_KEY_S) { if (n > 1) n--; changed = true; }
            else if (key == GLFW_KEY_E) { l++; changed = true; }
            else if (key == GLFW_KEY_D) { if (l > 0) l--; changed = true; }
            else if (key == GLFW_KEY_R) { m++; changed = true; }
            else if (key == GLFW_KEY_F) { m--; changed = true; }

            if (l >= n) l = n - 1;
            if (m > l) m = l;
            if (m < -l) m = -l;

            if (changed) {
                generateParticles(N);
                cout << "Orbital Updated -> n: " << n << " | l: " << l << " | m: " << m << "\n";
            }
        });
    }

    void render() {
        glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);

        mat4 proj = perspective(radians(45.0f), (float)w / (float)h, 0.1f, 500.0f);
        mat4 view = lookAt(camera.position(), camera.target, vec3(0, 1, 0));
        mat4 mvp = proj * view;

        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "uMVP"), 1, GL_FALSE, value_ptr(mvp));
        glUniform1f(glGetUniformLocation(shaderProgram, "uPointSize"), 6.0f);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, particles.size() * sizeof(Particle), particles.data());

        glBindVertexArray(VAO);
        glDrawArrays(GL_POINTS, 0, (GLsizei)particles.size());
        glBindVertexArray(0);
    }
};

//Main Loop
int main() {
    Engine engine;
    engine.setupCallbacks();

    generateParticles(N);

    cout << "===================================================\n";
    cout << " Quantum Orbital Probability Density Visualizer\n";
    cout << " Controls:\n";
    cout << "   [W/S] : Inc/Dec Principal Number (n)\n";
    cout << "   [E/D] : Inc/Dec Azimuthal Number (l)\n";
    cout << "   [R/F] : Inc/Dec Magnetic Number  (m)\n";
    cout << "   [Mouse Drag & Scroll] : Orbit / Zoom View\n";
    cout << "===================================================\n";

    float dt = 0.2f;
    double lastTime = glfwGetTime();

    while (!glfwWindowShouldClose(engine.window)) {
        double currentTime = glfwGetTime();
        float frameDt = (float)(currentTime - lastTime);
        lastTime = currentTime;

        // Kinematics: Update particle probability flow
        if (m != 0) {
            for (Particle& p : particles) {
                p.vel = calculateProbabilityFlow(p, m);
                p.pos += p.vel * dt;
            }
        }
        engine.render();
        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }
    glDeleteVertexArrays(1, &engine.VAO);
    glDeleteBuffers(1, &engine.VBO);
    glDeleteProgram(engine.shaderProgram);
    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}