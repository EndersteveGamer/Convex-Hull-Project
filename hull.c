#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct vec {
    double x;
    double y;
};

struct vecset {
    struct vec *data;
    size_t size;
    size_t capacity;
};

const int A = 2;

// Functions on vectors
double dot(struct vec *v1, struct vec *v2) {
    return v1->x * v2->x + v1->y * v2->y;
}

double cross(struct vec *v1, struct vec *v2, struct vec *v3) {
    return (v2->x - v1->x) * (v3->y - v1->y) - (v2->y - v1->y) * (v3->x - v1->x);
}

bool is_left_turn(struct vec *v1, struct vec *v2, struct vec *v3) {
    return cross(v1, v2, v3) < 0;
}

// Functions on vector sets
void vecset_create(struct vecset *self) {
    const int defaultCapacity = 10;
    self->capacity = defaultCapacity;
    self->size = 0;
    self->data = calloc(defaultCapacity, sizeof(struct vec));
}

void vecset_destroy(struct vecset *self) {
    free(self->data);
}

void vecset_grow(struct vecset *self) {
    self->capacity *= A;
    struct vec *new = calloc(self->capacity, sizeof(struct vec));
    memcpy(new, self->data, self->size * sizeof(struct vec));
    free(self->data);
    self->data = new;
}

int main() {
}
