#include <GL/glut.h>

#include "Renderer.h"

void initRendering() {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST); // Antialias the lines
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}


void render(
        int argc,
        char **argv,
        const char *title,
        int windowWidth,
        int windowHeight,
        void (*renderFn)(),
        void (*keyboardFn)(unsigned char, int, int)
) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(title);
    initRendering();

    glutDisplayFunc(renderFn);
    glutKeyboardFunc(keyboardFn);
    glutMainLoop();
}
