

#ifndef MOLECULAR_DYNAMICS_POC_RENDERER_H
#define MOLECULAR_DYNAMICS_POC_RENDERER_H


void render(
        int argc,
        char **argv,
        const char *title,
        int windowWidth,
        int windowHeight,
        void (*renderFn)(),
        void (*keyboardFn)(unsigned char, int, int)
);

#endif //MOLECULAR_DYNAMICS_POC_RENDERER_H
