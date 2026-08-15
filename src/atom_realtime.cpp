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

// Constants
const float a0 = 1.0f;
const double hbar = 1.0;
const double m_e = 1.0;
const double zmSpeed = 2.0;

// Quantum numbers & particle count (Optimized default for CPU laptops)
int n = 2, l = 1, m = 0;
int N_particles = 100000;

struct ParticleVertex {
    vec3 pos;
    vec4 color;
};

vector<ParticleVertex> particleData;
vector<vec3> particleVelocities;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<float> dis(0.0f, 1.0f);

// Color mapping (Heatmap)
vec4 heatmap_fire(float value) {
    value = std::clamp(value, 0.0f, 1.0f);
    const int num_stops = 6;
    const vec4 colors[num_stops] = {
        {0.0f, 0.0f, 0.0f, 1.0f},
        {0.5f, 0.0f, 0.99f, 1.0f},
        {0.8f, 0.0f, 0.0f, 1.0f},
        {1.0f, 0.5f, 0.0f, 1.0f},
        {1.0f, 1.0f, 0.0f, 1.0f},
        {1.0f, 1.0f, 1.0f, 1.0f}
    };

    float scaled_v = value * (num_stops - 1);
    int idx = static_cast<int>(scaled_v);
    int next_idx = std::min(idx + 1, num_stops - 1);
    float local_t = scaled_v - idx;

    return mix(colors[idx], colors[next_idx], local_t);
}

// Associated Laguerre polynomial evaluation
double evaluateLaguerre(int k, int alpha, double rho) {
    if (k <= 0) return 1.0;
    double Lm1 = 1.0 + alpha - rho;
    if (k == 1) return Lm1;
    
    double Lm2 = 1.0;
    double L = Lm1;
    for (int j = 2; j <= k; ++j) {
        L = ((2.0 * j - 1.0 + alpha - rho) * Lm1 - (j - 1.0 + alpha) * Lm2) / j;
        Lm2 = Lm1;
        Lm1 = L;
    }
    return L;
}

// Associated Legendre polynomial evaluation
double evaluateLegendre(int l_val, int m_val, double x) {
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
    if (l_val == abs_m) return Pmm;

    double Pm1m = x * (2 * abs_m + 1) * Pmm;
    if (l_val == abs_m + 1) return Pm1m;

    double Plm = Pm1m;
    for (int ll = abs_m + 2; ll <= l_val; ++ll) {
        Plm = ((2.0 * ll - 1.0) * x * Pm1m - (ll + abs_m - 1.0) * Pmm) / (ll - abs_m);
        Pmm = Pm1m;
        Pm1m = Plm;
    }
    return Plm;
}

// Dynamic CDF sampling for radial distance r
vector<double> r_cdf;
double rMax_cached = 0.0;
void rebuildRadialCDF(int cur_n, int cur_l) {
    const int steps = 2048;
    rMax_cached = 12.0 * cur_n * cur_n * a0;
    r_cdf.resize(steps);
    double dr = rMax_cached / (steps - 1);
    double sum = 0.0;

    int k = cur_n - cur_l - 1;
    int alpha = 2 * cur_l + 1;
    double norm = pow(2.0 / (cur_n * a0), 3) * tgamma(cur_n - cur_l) / (2.0 * cur_n * tgamma(cur_n + cur_l + 1));

    for (int i = 0; i < steps; ++i) {
        double r = i * dr;
        double rho = 2.0 * r / (cur_n * a0);
        double L = evaluateLaguerre(k, alpha, rho);
        double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, cur_l) * L;
        sum += r * r * R * R;
        r_cdf[i] = sum;
    }
    for (double& v : r_cdf) v /= sum;
}

double sampleR(mt19937& rng) {
    uniform_real_distribution<double> uniform_dis(0.0, 1.0);
    double u = uniform_dis(rng);
    int idx = lower_bound(r_cdf.begin(), r_cdf.end(), u) - r_cdf.begin();
    return idx * (rMax_cached / (r_cdf.size() - 1));
}

// Dynamic CDF sampling for theta
vector<double> theta_cdf;
void rebuildThetaCDF(int cur_l, int cur_m) {
    const int steps = 1024;
    theta_cdf.resize(steps);
    double dtheta = M_PI / (steps - 1);
    double sum = 0.0;

    for (int i = 0; i < steps; ++i) {
        double theta = i * dtheta;
        double x = cos(theta);
        double Plm = evaluateLegendre(cur_l, cur_m, x);
        sum += sin(theta) * Plm * Plm;
        theta_cdf[i] = sum;
    }
    for (double& v : theta_cdf) v /= sum;
}

double sampleTheta(mt19937& rng) {
    uniform_real_distribution<double> uniform_dis(0.0, 1.0);
    double u = uniform_dis(rng);
    int idx = lower_bound(theta_cdf.begin(), theta_cdf.end(), u) - theta_cdf.begin();
    return idx * (M_PI / (theta_cdf.size() - 1));
}

float samplePhi() {
    return 2.0f * (float)M_PI * dis(gen);
}

vec3 sphericalToCartesian(float r, float theta, float phi) {
    return vec3(r * sin(theta) * cos(phi), r * cos(theta), r * sin(theta) * sin(phi));
}

vec4 computeOrbitalColor(double r, double theta, int cur_n, int cur_l, int cur_m) {
    double rho = 2.0 * r / (cur_n * a0);
    double L = evaluateLaguerre(cur_n - cur_l - 1, 2 * cur_l + 1, rho);
    double norm = pow(2.0 / (cur_n * a0), 3) * tgamma(cur_n - cur_l) / (2.0 * cur_n * tgamma(cur_n + cur_l + 1));
    double R = sqrt(norm) * exp(-rho / 2.0) * pow(rho, cur_l) * L;

    double Plm = evaluateLegendre(cur_l, cur_m, cos(theta));
    double intensity = (R * R) * (Plm * Plm);

    return heatmap_fire((float)(intensity * 2.0 * pow(4.0, cur_n)));
}

vec3 calculateProbabilityFlow(const vec3& pos, int cur_m) {
    double r = length(pos);
    if (r < 1e-5) return vec3(0.0f);
    
    double theta = acos(std::clamp(pos.y / (float)r, -1.0f, 1.0f));
    double phi = atan2(pos.z, pos.x);

    double sinTheta = sin(theta);
    if (abs(sinTheta) < 1e-4) sinTheta = 1e-4;

    double v_mag = (hbar * cur_m) / (m_e * r * sinTheta);
    return vec3((float)(-v_mag * sin(phi)), 0.0f, (float)(v_mag * cos(phi)));
}

void generateParticles(int count) {
    rebuildRadialCDF(n, l);
    rebuildThetaCDF(l, m);

    particleData.resize(count);
    particleVelocities.resize(count);

    for (int i = 0; i < count; ++i) {
        float r = (float)sampleR(gen);
        float theta = (float)sampleTheta(gen);
        float phi = samplePhi();

        vec3 pos = sphericalToCartesian(r, theta, phi);
        vec4 col = computeOrbitalColor(r, theta, n, l, m);

        particleData[i] = { pos, col };
        particleVelocities[i] = vec3(0.0f);
    }
}

// Orbit Camera
struct Camera {
    vec3 target = vec3(0.0f);
    float radius = 40.0f;
    float azimuth = 0.0f;
    float elevation = (float)M_PI / 2.0f;
    float orbitSpeed = 0.005f;
    bool dragging = false;
    double lastX = 0.0, lastY = 0.0;

    vec3 position() const {
        float el = std::clamp(elevation, 0.01f, (float)M_PI - 0.01f);
        return vec3(radius * sin(el) * cos(azimuth), radius * cos(el), radius * sin(el) * sin(azimuth));
    }

    void processMouseMove(double x, double y) {
        if (dragging) {
            azimuth += (float)(x - lastX) * orbitSpeed;
            elevation -= (float)(y - lastY) * orbitSpeed;
            elevation = std::clamp(elevation, 0.01f, (float)M_PI - 0.01f);
        }
        lastX = x;
        lastY = y;
    }

    void processMouseButton(int button, int action, GLFWwindow* win) {
        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if (action == GLFW_PRESS) {
                dragging = true;
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if (action == GLFW_RELEASE) {
                dragging = false;
            }
        }
    }

    void processScroll(double yoffset) {
        radius -= (float)(yoffset * zmSpeed);
        if (radius < 2.0f) radius = 2.0f;
    }
} camera;

// Shader sources (Fast GL Points)
const char* particleVS = R"glsl(
    #version 330 core
    layout(location = 0) in vec3 aPos;
    layout(location = 1) in vec4 aColor;
    uniform mat4 view;
    uniform mat4 projection;
    out vec4 vColor;
    void main() {
        gl_Position = projection * view * vec4(aPos, 1.0);
        gl_PointSize = max(1.5, 80.0 / gl_Position.w);
        vColor = aColor;
    }
)glsl";

const char* particleFS = R"glsl(
    #version 330 core
    in vec4 vColor;
    out vec4 FragColor;
    void main() {
        vec2 coord = gl_PointCoord - vec2(0.5);
        if (dot(coord, coord) > 0.25) discard; // Smooth circular point sprites
        FragColor = vColor;
    }
)glsl";

const char* lineVS = R"glsl(
    #version 330 core
    layout(location = 0) in vec3 aPos;
    uniform mat4 view;
    uniform mat4 projection;
    void main() {
        gl_Position = projection * view * vec4(aPos, 1.0);
    }
)glsl";

const char* lineFS = R"glsl(
    #version 330 core
    uniform vec4 lineColor;
    out vec4 FragColor;
    void main() {
        FragColor = lineColor;
    }
)glsl";

GLuint compileShader(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    return s;
}

GLuint linkProgram(GLuint vs, GLuint fs) {
    GLuint p = glCreateProgram();
    glAttachShader(p, vs);
    glAttachShader(p, fs);
    glLinkProgram(p);
    return p;
}

int main() {
    if (!glfwInit()) return -1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(1024, 768, "Hydrogen Orbital & Probability Flow", NULL, NULL);
    if (!window) { glfwTerminate(); return -1; }
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) return -1;

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE); // Additive blending for glowing orbital look

    // Setup shaders
    GLuint pVS = compileShader(GL_VERTEX_SHADER, particleVS);
    GLuint pFS = compileShader(GL_FRAGMENT_SHADER, particleFS);
    GLuint particleProgram = linkProgram(pVS, pFS);

    GLuint lVS = compileShader(GL_VERTEX_SHADER, lineVS);
    GLuint lFS = compileShader(GL_FRAGMENT_SHADER, lineFS);
    GLuint lineProgram = linkProgram(lVS, lFS);

    GLint pViewLoc = glGetUniformLocation(particleProgram, "view");
    GLint pProjLoc = glGetUniformLocation(particleProgram, "projection");
    GLint lViewLoc = glGetUniformLocation(lineProgram, "view");
    GLint lProjLoc = glGetUniformLocation(lineProgram, "projection");
    GLint lColLoc  = glGetUniformLocation(lineProgram, "lineColor");

    // Grid geometry
    vector<float> gridVertices;
    float gridSize = 100.0f;
    int divisions = 20;
    float step = gridSize / divisions;
    for (int i = -divisions / 2; i <= divisions / 2; ++i) {
        gridVertices.insert(gridVertices.end(), { (float)i * step, 0.0f, -gridSize / 2.0f, (float)i * step, 0.0f, gridSize / 2.0f });
        gridVertices.insert(gridVertices.end(), { -gridSize / 2.0f, 0.0f, (float)i * step, gridSize / 2.0f, 0.0f, (float)i * step });
    }

    GLuint gridVAO, gridVBO;
    glGenVertexArrays(1, &gridVAO);
    glGenBuffers(1, &gridVBO);
    glBindVertexArray(gridVAO);
    glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
    glBufferData(GL_ARRAY_BUFFER, gridVertices.size() * sizeof(float), gridVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Initial Particle sampling
    generateParticles(N_particles);

    // Particle VAO / Dynamic VBO
    GLuint particleVAO, particleVBO;
    glGenVertexArrays(1, &particleVAO);
    glGenBuffers(1, &particleVBO);
    glBindVertexArray(particleVAO);
    glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
    glBufferData(GL_ARRAY_BUFFER, particleData.size() * sizeof(ParticleVertex), nullptr, GL_STREAM_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleVertex), (void*)offsetof(ParticleVertex, pos));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(ParticleVertex), (void*)offsetof(ParticleVertex, color));
    glEnableVertexAttribArray(1);

    // Callbacks
    glfwSetWindowUserPointer(window, &camera);
    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
        ((Camera*)glfwGetWindowUserPointer(win))->processMouseMove(x, y);
    });
    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int btn, int act, int) {
        ((Camera*)glfwGetWindowUserPointer(win))->processMouseButton(btn, act, win);
    });
    glfwSetScrollCallback(window, [](GLFWwindow* win, double, double yoff) {
        ((Camera*)glfwGetWindowUserPointer(win))->processScroll(yoff);
    });

    glfwSetKeyCallback(window, [](GLFWwindow*, int key, int, int action, int) {
        if (action != GLFW_PRESS) return;
        bool changed = false;
        if (key == GLFW_KEY_W) { n += 1; changed = true; }
        if (key == GLFW_KEY_S) { n = std::max(1, n - 1); changed = true; }
        if (key == GLFW_KEY_E) { l += 1; changed = true; }
        if (key == GLFW_KEY_D) { l = std::max(0, l - 1); changed = true; }
        if (key == GLFW_KEY_R) { m += 1; changed = true; }
        if (key == GLFW_KEY_F) { m -= 1; changed = true; }

        l = std::clamp(l, 0, n - 1);
        m = std::clamp(m, -l, l);

        if (changed) {
            cout << "State updated: n=" << n << " l=" << l << " m=" << m << "\n";
            generateParticles(N_particles);
        }
    });

    cout << "Controls:\n [W/S] n +/- | [E/D] l +/- | [R/F] m +/-\n Left-Drag: Orbit | Scroll: Zoom\n";

    float dt = 0.2f;
    while (!glfwWindowShouldClose(window)) {
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);

        glClearColor(0.02f, 0.02f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        mat4 projection = perspective(radians(45.0f), (float)width / (float)height, 0.1f, 1000.0f);
        mat4 view = lookAt(camera.position(), camera.target, vec3(0, 1, 0));

        // 1. Draw Grid
        glUseProgram(lineProgram);
        glUniformMatrix4fv(lViewLoc, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(lProjLoc, 1, GL_FALSE, value_ptr(projection));
        glUniform4f(lColLoc, 0.2f, 0.25f, 0.35f, 0.4f);
        glBindVertexArray(gridVAO);
        glDrawArrays(GL_LINES, 0, (GLsizei)gridVertices.size() / 3);

        // 2. Update probability flow physics (phi rotation)
        if (m != 0) {
            for (size_t i = 0; i < particleData.size(); ++i) {
                vec3& p = particleData[i].pos;
                double r = length(p);
                if (r > 1e-4) {
                    double theta = acos(std::clamp(p.y / (float)r, -1.0f, 1.0f));
                    vec3 vel = calculateProbabilityFlow(p, m);
                    vec3 temp_pos = p + vel * dt;
                    double new_phi = atan2(temp_pos.z, temp_pos.x);
                    p = sphericalToCartesian((float)r, (float)theta, (float)new_phi);
                }
            }
        }

        // 3. Render Particles in a SINGLE draw call
        glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
        glBufferData(GL_ARRAY_BUFFER, particleData.size() * sizeof(ParticleVertex), particleData.data(), GL_STREAM_DRAW);

        glUseProgram(particleProgram);
        glUniformMatrix4fv(pViewLoc, 1, GL_FALSE, value_ptr(view));
        glUniformMatrix4fv(pProjLoc, 1, GL_FALSE, value_ptr(projection));

        glBindVertexArray(particleVAO);
        glDrawArrays(GL_POINTS, 0, (GLsizei)particleData.size());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}